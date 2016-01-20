##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=libl2cbitstream
ConfigurationName      :=Debug
WorkspacePath          := "/home/pasi/sandbox/pmiettinen/l2cbitstream"
ProjectPath            := "/home/pasi/sandbox/pmiettinen/l2cbitstream/libl2cbitstream"
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Pasi Miettinen
Date                   :=20/01/16
CodeLitePath           :="/home/pasi/.codelite"
LinkerName             :=/usr/bin/g++
SharedObjectLinkerName :=/usr/bin/g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=./$(ProjectName).so
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="libl2cbitstream.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)/usr/include/python2.7 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)swiftnav $(LibrarySwitch)python2.7 
ArLibs                 :=  "libswiftnav" "python2.7" 
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -g -fPIC $(Preprocessors)
CFLAGS   :=  -g -fPIC $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/l2cbitstream.c$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(SharedObjectLinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)
	@$(MakeDirCommand) "/home/pasi/sandbox/pmiettinen/l2cbitstream/.build-debug"
	@echo rebuilt > "/home/pasi/sandbox/pmiettinen/l2cbitstream/.build-debug/libl2cbitstream"

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/l2cbitstream.c$(ObjectSuffix): l2cbitstream.c $(IntermediateDirectory)/l2cbitstream.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/pasi/sandbox/pmiettinen/l2cbitstream/libl2cbitstream/l2cbitstream.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/l2cbitstream.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/l2cbitstream.c$(DependSuffix): l2cbitstream.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/l2cbitstream.c$(ObjectSuffix) -MF$(IntermediateDirectory)/l2cbitstream.c$(DependSuffix) -MM "l2cbitstream.c"

$(IntermediateDirectory)/l2cbitstream.c$(PreprocessSuffix): l2cbitstream.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/l2cbitstream.c$(PreprocessSuffix) "l2cbitstream.c"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


