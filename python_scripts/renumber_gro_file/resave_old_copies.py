import os
import re

def store_old_copy_gromacs_style(old_file, search_dir='.', look_for_extension=True): #look_for_extension=True means that the script looks for an extension code at the end of the filename   

    #could insert here code to find path of old_file (if it isn't in cwd), so that search_dir will be automatically located
 
    dirlist = os.listdir(search_dir)
    if old_file not in dirlist:
        return None
    print 'Old copy of {} located'.format(old_file)
    
    #find out whether filename has format <file>.XXX  (eg. <file>.gro) or just <file>, and set up variable 'old_file_stem', which, will be <file>, most likely, and old_file_extension
    if look_for_extension==True:
        if re.search(r'(.+)\.\w{3}', old_file):
            extension = True
            old_file_stem = re.sub(r'(.+)\.\w{3}', r'\1', old_file)
            old_file_extension = re.sub(r'.+(\.\w{3})', r'\1', old_file)
        else:
            extension = False
            old_file_stem = filename
            old_file_extension = ''
    else:
        old_file_stem = filename
        old_file_extension = ''
     
    # find correct version number        
    last_version = 0
    for filename in dirlist:
        m = re.match(r'#%s\.(\d+)%s#' % (old_file_stem, old_file_extension), filename)
        if m:
            version = int(m.group(1))
            if version > last_version:
                last_version = version     
    
    version = last_version + 1
    old_file_newname = '#%s.%d%s#' % (old_file_stem, version, old_file_extension)
    os.rename(old_file, old_file_newname)
    print 'Old copy of %s has been backed up to %s' % (old_file, old_file_newname)



