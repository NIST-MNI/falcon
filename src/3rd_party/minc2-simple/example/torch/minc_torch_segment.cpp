#include <torch/torch.h>
#include <torch/script.h> // One-stop header


#include "minc2-simple.h"

#include <iostream>
#include <memory>

// for get_opt_long
#include <getopt.h>

void show_usage(const char *name)
{
  std::cerr 
    << "Usage: "<<name<<" <script.pth> <input> <output> [debug]" << std::endl
    << "Optional parameters:" << std::endl
    << "\t--channels <n> add more blank input channels, set them to 38.81240207 if missing for now" << std::endl
    << "\t--add <n> add more input channels " << std::endl
    << "\t--mask <file> use binary mask to restrict application of the network (TODO)" << std::endl
    << "\t--verbose be verbose" << std::endl
    << "\t--clobber clobber the output files" << std::endl
    << "\t--patch <n> patch size" << std::endl
    << "\t--stride <n> stride, should be less then patch" << std::endl
    << "\t--crop <n>   crop n voxels from the edge, should be less then patch/2" << std::endl;
}

torch::Tensor load_float_tensor(minc2_file_handle h)
{
    struct minc2_dimension * _dims;
    minc2_get_representation_dimensions(h, &_dims);

    std::cout<<"load_float_tensor:" << _dims[0].length << "," << _dims[1].length << "," << _dims[2].length << std::endl;

    torch::Tensor output  = torch::empty( { _dims[2].length, _dims[1].length, _dims[0].length }, torch::kFloat32);
    if(minc2_load_complete_volume(h, output.data_ptr(), MINC2_FLOAT)!=MINC2_SUCCESS)
    {
            std::cerr << "Error reading data from handle";
            throw std::exception();
    }
    return output;
}

template<class tensor_type> void save_tensor(const char* out_file,torch::Tensor t, minc2_file_handle ref,int minc_data_type)
{
    minc2_file_handle o = minc2_allocate0();
    struct minc2_dimension * store_dims;
    minc2_get_store_dimensions(ref, &store_dims);
    std::cout<<"save_tensor: store: " << store_dims[0].length << ","<< store_dims[1].length << "," << store_dims[2].length << std::endl;

    minc2_define(o, store_dims, minc_data_type, minc_data_type);

    if(minc2_create(o, out_file)!=MINC2_SUCCESS)
    {
        std::cerr << "Can't open " << out_file << "for writing" << std::endl;
        throw std::exception();
    }
    minc2_setup_standard_order(o);

    if(minc2_save_complete_volume(o, t.data_ptr<tensor_type>(), minc_data_type)!=MINC2_SUCCESS)
    {
        std::cerr << "Error writing data to " << out_file << std::endl;
        throw std::exception();
    }
    minc2_close(o);
    minc2_free(o);
}


int main(int argc,char **argv) 
{
    int clobber=0;
    int verbose=0;
    std::string mask_f;
    int channels=1;
    int c;
    int patch_sz = 80;
    int crop = 4;
    int stride = (patch_sz-crop*2)/2;
    std::vector<std::string> add_f;

    // read the arguments
    static struct option long_options[] =
    {
        {"verbose",   no_argument, &verbose, 1},
        {"quiet",     no_argument, &verbose, 0},
        {"clobber",   no_argument, &clobber, 1},
        {"mask",      required_argument, 0, 'm'},
        {"channels",  required_argument, 0, 'c'},
        {"patch",     required_argument, 0, 'p'},
        {"stride",    required_argument, 0, 's'},
        {"crop",      required_argument, 0, 'R'},
        {"add",       required_argument, 0, 'A'},
        {0, 0, 0, 0}
    };

    for (;;)
    {
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "v", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 0:
            break;
        case 'c':
            channels=atoi(optarg);
            break;
        case 'p':
            patch_sz=atoi(optarg);
            break;
        case 's':
            stride=atoi(optarg);
            break;
        case 'R':
            crop=atoi(optarg);
            break;
        case 'm':
            mask_f=optarg;
            break;
        case 'A':
            add_f.push_back(optarg);
            break;
        case '?':
                /* getopt_long already printed an error message. */
        default:
                show_usage(argv[0]);
                return 1;
        }
    }
    if ((argc - optind) < 3)
    {
        show_usage(argv[0]);
        return 1;
    }
    bool debug_output = (argc - optind)>3;
    
    try {
        int ndim;
        int nelement;
        int patch_sz_eff = patch_sz - crop*2;
        struct minc2_dimension * _dims;
        torch::NoGradGuard no_grad; // disable gradient calculation (?)

        // Deserialize the ScriptModule from a file using torch::jit::load().
        torch::jit::script::Module module = torch::jit::load(argv[optind]);
        module.eval();

        // dump module
        // module.dump(true,true,false);

        // load minc file
        minc2_file_handle h = minc2_allocate0();

        if(minc2_open(h, argv[optind+1])!=MINC2_SUCCESS)
        {
            std::cerr << "Can't open " << argv[optind]  << " for reading" << std::endl;
            throw std::exception();
        }

        minc2_setup_standard_order(h);
        minc2_ndim(h, &ndim);
        minc2_nelement(h, &nelement);
        minc2_get_representation_dimensions(h, &_dims);

        torch::Tensor input = load_float_tensor(h);

        torch::Tensor output;
        torch::Tensor output_d;

        torch::Tensor ones;
        torch::Tensor weights;

        // convert to 5D format
        input = input.unsqueeze(0).unsqueeze(0);

        // add more input channels 
        // TODO: read them from file
        if(channels>1) {
            torch::Tensor dummy = torch::full( { _dims[2].length, _dims[1].length, _dims[0].length }, 38.81240207);
            dummy = dummy.unsqueeze(0).unsqueeze(0);
            std::vector<torch::Tensor> _inputs={input};
            for(int i=1;i<channels;i++)
                _inputs.push_back(dummy);
            // append new channels
            input = torch::cat(_inputs,1);
        } else if(!add_f.empty()) {
            std::vector<torch::Tensor> _inputs={input};
            for(auto a:add_f)
            {
                minc2_file_handle hh = minc2_allocate0();
                if(minc2_open(hh, a.c_str())!=MINC2_SUCCESS)
                {
                    std::cerr << "Can't open " << a << " for reading" << std::endl;
                    throw std::exception();
                }
                minc2_setup_standard_order(hh);

                _inputs.push_back(load_float_tensor(hh).unsqueeze(0).unsqueeze(0));
                minc2_close(hh);
                minc2_free(hh);
            }
            input = torch::cat(_inputs,1);
        }

        // TODO: queue patches in the batch dimension?
        for(int z=0; z< ceil((double)_dims[2].length / stride ); z ++)
        {
            for(int y=0; y< ceil((double)_dims[1].length / stride ); y ++)
            {
                for(int x=0; x< ceil((double)_dims[0].length / stride ); x ++)
                {
                    torch::NoGradGuard no_grad;
                    int _x = x*stride;
                    int _y = y*stride;
                    int _z = z*stride;

                    _x = ( _x + patch_sz) > _dims[0].length ? _dims[0].length - patch_sz:_x;
                    _y = ( _y + patch_sz) > _dims[1].length ? _dims[1].length - patch_sz:_y;
                    _z = ( _z + patch_sz) > _dims[2].length ? _dims[2].length - patch_sz:_z;

                    std::cout << _x << "," << _y << "," << _z <<std::endl;
                    // Create a vector of inputs.
                    std::vector<torch::jit::IValue> inputs;

                    auto inp = input.slice(2,_z,_z+patch_sz).slice(3,_y,_y+patch_sz).slice(4, _x, _x+patch_sz);
                    // extract patch
                    inputs.push_back( inp );

                    // Execute the model and turn its output into a tensor.
                    auto out = torch::log_softmax( module.forward(inputs).toTensor(), /*dim=*/1 );

                    // drop batch dimension
                    out = out.squeeze(0);

                    if(!output.defined()) // neeed to allocate
                    {
                        std::cout<<"out:"<< out.size(0)<<","<< out.size(1)<<","<< out.size(2)<<","<< out.size(3) << std::endl;
                        std::cout<<"Allocating output:" << out.size(0) << ","<< _dims[0].length << ","<< _dims[1].length << "," << _dims[2].length << std::endl;
                        std::cout<<"Allocating ones:"   << out.size(0) << ","<< patch_sz << "," << patch_sz << "," << patch_sz << std::endl;

                        output  = torch::zeros( { out.size(0), _dims[2].length, _dims[1].length, _dims[0].length });
                        weights = torch::zeros( { out.size(0), _dims[2].length, _dims[1].length, _dims[0].length });
                        ones    = torch::ones(  { out.size(0), patch_sz_eff,       patch_sz_eff, patch_sz_eff    });
                    }

                    // remove boundary voxels
                    output.slice(1, _z+crop ,_z+crop+patch_sz_eff).slice(2,_y+crop,_y+crop+patch_sz_eff).slice(3,_x+crop,_x+crop+patch_sz_eff) += 
                        out.slice(1,crop,crop+patch_sz_eff).slice(2,crop,crop+patch_sz_eff).slice(3,crop,crop+patch_sz_eff);
                    // account for overlapping samples
                    weights.slice(1,_z+crop ,_z+crop+patch_sz_eff).slice(2,_y+crop,_y+crop+patch_sz_eff).slice(3,_x+crop,_x+crop+patch_sz_eff) += ones;

                }
            }
        }
        // HACK to avoid devide by zero on the border
        torch::Tensor undefined_roi = weights<1.0;
        weights.masked_fill_(undefined_roi,1.0);

        output= output / weights; // normalize
        output.exp_(); // apply exponent

        if(debug_output) // output 1st class , usually BG
            output_d = output.slice(0,0,1).squeeze(0);

        // extract prediction
        auto out_ = output.argmax(0).to(torch::kU8); 
        out_.masked_fill_(undefined_roi.slice(0,0,1).squeeze(0),0); // assignd default label to voxels with undefined results

        save_tensor<u_char>(argv[optind+2],out_,h,MINC2_UBYTE);

        if(debug_output)
        {
            save_tensor<float>(argv[optind+3],output_d,h,MINC2_FLOAT);
        }

        minc2_close(h);

        /*deallocate*/
        minc2_free(h);
    }
    catch (const c10::Error& e) 
    {
        std::cerr << "error:"<< e.what() << std::endl;
        return -1;
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Error " << std::endl;
        return -1;
    }
    std::cout << "ok\n";
    return 0;
}
