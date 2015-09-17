# set up build directory
builddir=${1:-build/} # by default use build/
echo "Build directory: ${builddir}"
mkdir -p $builddir
( cd $builddir && cmake -DCMAKE_BUILD_TYPE=Coverage .. )

hooks="`git rev-parse --git-dir`/hooks"
echo "Installing hook at: ${hooks}"
# insert value of builddir as first line of pre-commit
echo "subdir=$builddir" | cat - checks/pre-commit > "$hooks/pre-commit"
