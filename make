#!/bin/bash

#	Name: make.sh
#	Author: Johan Hidding
#	Licence: do whatever you like with this, there's no waranty.
# -----------------
#	Description: this simple script should be able to compile and
#	link a simple but well structured c/c++ project.
# -----------------
#	Details: compiles every .cc file found in all subdirectories. 
#	the script checks the modification times of the source files
#	and their header dependencies against the target files.
#	finally all .o files are linked to $target.
#	if you want to exclude files from compiling, put them in a
#	hidden directory (or somewhere outside the directory
#	structure, which is even better)

target="dmt"
objdir=".obj"
LDFLAGS="-lm -lfftw3 -fopenmp -lgsl -lgslcblas"
CFLAGS="-O2 -Wall -g -std=c++0x -fopenmp"

CC="g++"
ext=".cc"
ECHO="echo -e"

#/>==--------                                  --------==>\#
#>=----- you shouldn't need to edit below this line -----=<#
#\>==--------                                  --------==>/#

DIRS=`find . -maxdepth 2 -type d -name '[^\.]*'`
CCFILES=`find . -maxdepth 3 -wholename "./[^.]*$ext"`

case "$TERM" in
	dumb)
		prettyprint() {
			$ECHO " * $3 "
			if $1; then
				# $ECHO "\t[done]"
				return 1
			else
				$ECHO "\t[failed]"
				return 0
			fi
		} ;;
	*)
		prettyprint() {
			$ECHO "\033[$2m*\033[m $3"
			if $1; then
				$ECHO "\r\033[A\033[50C[\033[32mdone\033[m]"
				return 1
			else
				$ECHO "[\033[31mfailed\033[m]"
				return 0
			fi
		} ;;
esac

checknewer() {
	local k
	if [ ! -e $1 ]; then
		return 0
	fi

	for k in $2; do
		if [ $k -nt $1 ]; then
			return 0
		fi
	done
	return 1
}

compile() {
	objf=$objdir/$(basename $1 $ext).o
	deps=`$CC -MM $1 $CFLAGS | sed -e '{ s/^.*: //; s/\\\//; s/^ *// }'`
	if checknewer $objf "$deps make"; then
		if prettyprint "$CC -c $CFLAGS $1 -o $objf" 34 "Compiling $1 ... "; then
			exit 1
		fi
	fi
}

if [ ! -e $objdir ]; then
	mkdir $objdir
fi

case "$1" in
	single)
		compile $2 ;;
		
	all)
		for f in $CCFILES; do
			compile $f
		done

		if checknewer $target "`ls $objdir/*.o`"; then
			if prettyprint "$CC $objdir/*.o -o $target $LDFLAGS" 36 "Linking ..."; then
				exit 1
			fi
		else
			echo "$target allready is up to date."
		fi ;;

	run)
		./make.sh -all
		exec ./$target ;;

	clean)
		rm -rf $target $objdir
		find . -name '*~' -exec rm {} \; ;;

	*)
		echo "This is a make script, written in bash. It is only good"
		echo "for compiling well-structured C++ programs. Read the"
		echo "source for help on configuring the script."
		echo
		echo "run> ./make [command]"
		echo "where [command] ∈ {single, all, run, clean}."
		echo
		echo "'make single' expects one more argument giving a .cc file"
		echo "that you want to compile."
esac

