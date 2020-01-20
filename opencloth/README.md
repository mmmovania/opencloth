Open Cloth solution for C++ projects, except for the CUDA project.

The C++ projects are included in this one solution, so it's easy to compile all projects at once.

A common properties sheet has been added to the C++ projects, so changing paths to the include directories or to the library directories can be done in one spot and will apply to all projects.

Open Cloth C++ projects depend on freeglut, glew, and glm, and those dependencies are located in the 'dep' folder.  The common property sheet has relative paths to the dep folder, so no changes should be needed to compile the solution.

Clone or download the github repository, open 'opencloth.sln' in Visual Studio, select either 'debug' or 'release' targets, and click on 'rebuild solution'.

The included Visual Studio C++ projects were upgraded to VC++ 2010 Express Edition from VC++ 2008.  VC++ 2010 projects can be opened with Visual Studio 2012, 2015, 2017, or 2019.  Upon first opening, newer versions of Visual Studio will perform a one-way upgrade on the projects.  Then proceed with 'rebuild solution' as above.

