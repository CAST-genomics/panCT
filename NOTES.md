1. How do we want to include the gbz-base source code in the repo? Options: git submodule or install as setup step
2. Do we want to automatically skip compilation of the cython module if it doesn't work? Or should we just fail/error explicitly. The former could be helpful if we want to implement some kind of slower alternative
