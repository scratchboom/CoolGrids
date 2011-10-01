#pragma once

#include <string>
#include <string.h>//for memcpy
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include <set>

#include <sys/time.h>
#include <stdint.h>

#include <omp.h>
#include <mpi.h>




#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>

#include "IO.hpp"
#include "MpiWrapper.hpp"

#define cimg_display 0
#include <CImg.h>

#include "DefineSwitchKeys.hpp"

#include "readonly.hpp"


#include "Memory.hpp"
#include "Utils.hpp"
#include "Functions.hpp"
#include "StringFunctions.hpp"
#include "Solving.hpp"

#include "Timer.hpp"

#include "ByteBuffer.hpp"

#include "BmpImage.hpp"

#include "Grid1D.hpp"
#include "Grid2D.hpp"
#include "Grid3D.hpp"

#include "TimedGrid1D.hpp"
#include "TimedGrid2D.hpp"
#include "TimedGrid3D.hpp"

#include "ColorMap.hpp"

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

#include "GnuPlotSaver1D.hpp"
#include "GnuPlotSaver2D.hpp"
#include "CsvSaver.hpp"
#include "CsvTableSaver.hpp"
#include "BmpSaver.hpp"
#include "VtiSaver.hpp"
#include "CImgSaver2D.hpp"

#include "Units.hpp"

#include "ElectromagneticField.hpp"

#include "RiemannSolver.hpp"
#include "Hydrodynamics.hpp"

#include "HdVectors.hpp"

