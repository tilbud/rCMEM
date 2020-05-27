# Script to check for dependent functions in the JimH-dev scripts.
# Prints the function name followed by lines in which the function is called in the other scripts.
# functionnames.txt is a space delimited text file containing the function names.

for functionname in $(< functionnames.txt)
do
	echo $functionname >> Jim_Scripts_dependencies.txt
	grep $functionname *.R >> Jim_Scripts_dependencies.txt
	echo "--------------" >> Jim_Scripts_dependencies.txt
	echo "       " >> Jim_Scripts_dependencies.txt
done
