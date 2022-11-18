/* This is a hack*/
/* To be able to use standard libnifti, when it's bundled with libminc*/

int nifti_io_problem_code_kunio_dummy = 0;

/* for calling from some main program */
/* WARNING: this will not do anything usefull!*/

int nifti_io_get_problem_code_kunio() {
  return nifti_io_problem_code_kunio_dummy;
}

void nifti_io_set_problem_code_kunio(int i) {
  nifti_io_problem_code_kunio_dummy = i;
}
