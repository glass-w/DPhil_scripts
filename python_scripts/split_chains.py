from pymol import cmd, CmdException

def split_chains(selection='(all)', prefix=None):
    '''
DESCRIPTION

    Create a single object for each chain in selection

SEE ALSO

    split_states, http://pymolwiki.org/index.php/Split_object
    '''
    count = 0
    models = cmd.get_object_list('(' + selection + ')')
    for model in models:
        for chain in cmd.get_chains('(%s) and model %s' % (selection, model)):
            if chain == '':
                chain = "''"
            count += 1
            if not prefix:
                name = '%s_%s' % (model, chain)
            else:
                name = '%s%04d' % (prefix, count)
            cmd.create(name, '(%s) and model %s and chain %s' % (selection, model, chain))
        cmd.disable(model)

cmd.extend('split_chains', split_chains)

