#pragma once




#ifdef GRIDS_DONT_CHECK_INDICES
    #pragma message( "Grids++ : index checking is OFF" )
#else
    #define GRIDS_CHECK_INDICES
    #pragma message( "Grids++ : index checking is ON" )
#endif






//#undef ENABLE_OMP
/*
#ifdef ENABLE_OMP
    #define PARALLEL_FOR                  parallel for
    #define PARALLEL_FOR_COLLAPSE(n)      parallel for collapse(n)
#else
    #define PARALLEL_FOR                  nothing
    #define PARALLEL_FOR_COLLAPSE(n)      nothing
#endif
*/
