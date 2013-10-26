/// matlab like timer

// start the timer
void tic();

// return time since last call to tic()
double toc0();

// return time since last call to tic()
// it's show
//   `you string`: t = %.6f s.\n
double toc(const char *st);

// return time since last call to tic()
// if call as toc(), show
//   Elapsed time is %.6f seconds.\n
// if call as toc("t = %.3f\n"), show
//   t = %.3f\n
double tocs(const char *st);
