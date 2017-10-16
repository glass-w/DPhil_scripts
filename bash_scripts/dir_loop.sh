ORIG=`pwd`
i=0

for d in ./*/ ; do

        cd "$d"
	mkdir input output analysis

        BASHSCRIPTS=/biggin/b123/sedm5059/SCRIPTS/bash_scripts
	PYTHONSCRIPTS=/biggin/b123/sedm5059/SCRIPTS/python_scripts
	VMDSCRIPTS=/biggin/b123/sedm5059/SCRIPTS/vmd_scripts

	

	((i++))

	
	mkdir input output         
         
        cd $ORIG

done

