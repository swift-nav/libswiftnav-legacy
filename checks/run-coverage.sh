# requires gcovr and diff-cover
# pip install gcovr diff-cover
set -e
gcovr -r .. --object-directory . --xml -o ../.coverage.xml
gcovr -r .. --object-directory . --html -o coverage.html
( cd .. && diff-cover .coverage.xml --fail-under=85 --html-report .diffreport.html )
mv ../.diffreport.html ./diffreport.html
