nvcc -arch=sm_60 SolverCLI.cu -I../include -I"C:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\include" --std c++17 -L"C:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\lib\x64_windows_msvc14\stat_mda" cplex2010.lib -Xcompiler="/permissive- /GS /GL /W3 /Gy /Zc:wchar_t /Zi /Gm- /O2 /sdl /Zc:inline /fp:precise /D "NDEBUG" /D "_CONSOLE" /D "_UNICODE" /D "UNICODE" /WX- /Zc:forScope /Gd /Oi /MD /openmp /std:c++17 /FC"



