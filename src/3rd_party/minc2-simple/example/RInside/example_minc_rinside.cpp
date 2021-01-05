#include <RInside.h>                    // for the embedded R via RInside
#include <unistd.h>
#include <getopt.h>
#include "minc2-simple.h"

void show_usage (const char * prog)
{
  std::cerr<<"Usage:"<<prog<<" in.csv out_base "<<std::endl
      <<"[ "<<std::endl
      <<"  --clobber clobber output file(s)" <<std::endl
      <<"  --mask <f> use mask" <<std::endl
      <<"]"<<std::endl; 
}

int main(int argc, char *argv[]) {
    int verbose=0,clobber=0;
    int order=2;
    std::string mask_f;
    
    static struct option long_options[] = {
        {"verbose", no_argument,       &verbose, 1},
        {"quiet",   no_argument,       &verbose, 0},
        {"clobber", no_argument,       &clobber, 1},
        {"mask", required_argument,     0, 'm'},
        {0, 0, 0, 0}
    };
    
    for (;;) {
        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long (argc, argv, "vo:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1) break;

        switch (c)
        {
        case 0:
            break;
        case 'm':
            mask_f=optarg;break;
        case 'v':
            std::cout << "Version: 0.1" << std::endl;
            return 0;
        case '?':
            /* getopt_long already printed an error message. */
        default:
            show_usage (argv[0]);
            return 1;
        }
    }

    if((argc - optind) < 2) {
        show_usage (argv[0]);
        return 1;
    }
    std::string input_csv=argv[optind];
    std::string output_base=argv[optind+1];
    
    try {
        // create an embedded R instance
        RInside R(argc, argv,true,false,false ); //Verbose
        R["input"]=input_csv;

        //std::string r_code = "myRlist[[1]] = 42; myRlist[[2]] = 42.0; myRlist[[3]][2] = 42; myRlist";

        Rcpp::DataFrame csv_contents = R.parseEval("read.csv(input,header=F)");
        
        std::cout<<"Read csv file with :"<<csv_contents.size()<<" columns"<<std::endl;
        
        Rcpp::CharacterVector in_files=Rcpp::as<Rcpp::CharacterVector>(csv_contents[0]);

        
        std::cout<<"rows:"<<csv_contents.nrows()<<std::endl;
        /*opening minc files*/
        std::vector<minc2_file_handle> handles;
        for(Rcpp::CharacterVector::iterator it=in_files.begin();it!=in_files.end();++it)
        {
            std::cout<<"File:"<<*it<<std::endl;
            minc2_file_handle h=minc2_allocate0();
            if(minc2_open(h,*it)!=MINC2_SUCCESS)
                throw "Can't open file";
            handles.push_back(h);
        }
        
        
        R["voxel_data"]=csv_contents;
        R.parseEval("voxel_data[1]<-NULL");
        R.parseEval("voxel_data$voxel<-0.0");
        Rcpp::DataFrame voxel_data=R["voxel_data"];
        
        
        std::vector<minc2_file_handle> out_handles;
        minc2_file_handle o=minc2_allocate0();
        struct minc2_dimension *store_dims;
        minc2_get_store_dimensions(handles[0],&store_dims);
        minc2_define(o,store_dims,MINC2_DOUBLE, MINC2_DOUBLE);
        std::string out_fstat=output_base+"_fstat.mnc";
        minc2_create(o,out_fstat.c_str());
        out_handles.push_back(o);
        
        minc2_file_iterator_handle input_minc_it=minc2_iterator_allocate0();
        minc2_file_iterator_handle output_minc_it=minc2_iterator_allocate0();
        
        minc2_multi_iterator_input_start(input_minc_it,&handles[0],MINC2_DOUBLE,handles.size());
        minc2_multi_iterator_output_start(output_minc_it,&out_handles[0],MINC2_DOUBLE,out_handles.size());
        
        //Rcpp::NumericVector voxel=voxel_data["voxel"];
        //voxel[1]=1.0;
        //voxel_data["voxel"]=voxel;
        
        //R.parseEval("print(voxel_data)");
        //iterate over input files
        std::vector<double> voxels(handles.size());
        int cnt=0;
        do {
            cnt++;
            // read voxels into array
            minc2_iterator_get_values(input_minc_it,&voxels[0]);
            
            voxel_data["voxel"]=Rcpp::NumericVector(voxels.begin(),voxels.end());
            R.parseEval("res=summary(lm(voxel~.,data=voxel_data))");
            Rcpp::NumericVector qq=R.parseEval("res$fstatistic");
            double val=qq[0];
            //double val=Rcpp::as<double>(R.parseEval("mean(voxel_data$voxel)"));
            
            minc2_iterator_put_values(output_minc_it,&val);
            
            minc2_iterator_next(output_minc_it);
            if((cnt%1000)==1)
                std::cout<<cnt<<"\t"<<std::flush;
            
        } while(minc2_iterator_next(input_minc_it)==MINC2_SUCCESS);
        std::cout<<"Done:"<<cnt<<std::endl;
        
        minc2_iterator_free(output_minc_it);
        minc2_iterator_free(input_minc_it);
        for(std::vector<minc2_file_handle>::iterator it=handles.begin();it!=handles.end();++it)
        {
            minc2_close(*it);
            minc2_free(*it);
        }
        for(std::vector<minc2_file_handle>::iterator it=out_handles.begin();it!=out_handles.end();++it)
        {
            minc2_close(*it);
            minc2_free(*it);
        }
    } catch(std::exception& ex) {
        std::cerr << "Exception caught: " << ex.what() << std::endl;
    } catch(...) {
        std::cerr << "Unknown exception caught" << std::endl;
    }

    
}
