
#--------------------------------------------------------------------------
# Set Trilinos installation directory
#--------------------------------------------------------------------------
TRILINOS_INSTALL_DIR = ${HOME}/Code/BUILDS/OPT

#--------------------------------------------------------------------------
# If needed, set extra include and lib macros
#--------------------------------------------------------------------------
CLIENT_EXTRA_INCLUDES = 
CLIENT_EXTRA_LIBS = 

trackP.exe: trackP.o $(EXTRA_OBJS)

EXTRA_OBJS = RadiationModel.o StoichModel.o DiscretizationControl.o TransportModel.o OutputManager.o

RadiationModel.o: RadiationModel.cpp RadiationModel.hpp

StoichModel.o: StoichModel.cpp StoichModel.hpp

TransportModel.o: TransportModel.cpp TransportModel.hpp

DiscretizationControl.o: DiscretizationControl.cpp DiscretizationControl.hpp

OutputManager.o: OutputManager.cpp OutputManager.hpp

trackP.o: trackP.cpp






#----------------------------------------------------------------------------------------------
# THE REMAINING LINES WILL NOT NORMALLY NEED TO BE CHANGED. 
#
# Cases where you'll need to edit them include 
# (1) Your application needs some compiler flags not used in the Trilinos
#     build. Be careful to avoid inconsistency: for instance, mixing C++ 
#     libraries built with and without STL checking can cause segfaults.
# (2) Your application must link to some 3rd party libraries not among those 
#     specified when building Trilinos. You'll need to add them to the
#     linker command line. You may need to add their locations
#     to the library search path, and to the rpath if you're using shared
#     libraries.  
#-----------------------------------------------------------------------------------------------

# Include the Trilinos export makefile from package=Sundance.
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.Sundance

# Add the Trilinos installation directory to the search paths
# for libraries and headers
LIB_PATH = $(TRILINOS_INSTALL_DIR)/lib

INCLUDE_PATH = $(TRILINOS_INSTALL_DIR)/include $(CLIENT_EXTRA_INCLUDES)

# Set the C++ compiler and flags to those specified in the export makefile
CXX = $(Sundance_CXX_COMPILER)

CXXFLAGS = $(Sundance_CXX_FLAGS)

# Add the Trilinos libraries, search path, and rpath to the 
# linker command line arguments 
LIBS = $(CLIENT_EXTRA_LIBS) $(Sundance_SHARED_LIB_RPATH_COMMAND) \
  $(Sundance_LIBRARIES) \
  $(Sundance_TPL_LIBRARIES) $(Sundance_EXTRA_LD_FLAGS) 

# Rules for building executables and objects. 
%.exe : %.o $(EXTRA_OBJS)
	$(CXX) -o $@ $(LDFLAGS) $(CXXFLAGS) $< $(EXTRA_OBJS) -L$(LIB_PATH) $(LIBS)

%.o : %.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) -I$(INCLUDE_PATH) $(Sundance_TPL_INCLUDES) $<


