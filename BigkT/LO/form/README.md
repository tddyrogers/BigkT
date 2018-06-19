# FORM

## quick guide
- Run the script: ``install``. Some comments are available there.
- Copy the output of qgraf, e.g. the data folder here.
- The modify the ```template.frm``` which will be used by form 
- Run the ``form.py``

## What does form.py do? 
- It will take the ``qgraf.obj``` file and construct a form code
  using the ```template.frm``` to perform traces, replacements etc.
- The results will be placed at data.
  + ```form-mode-0.obj```: nothing will be done 
  + ```form-mode-1.obj```: traces won't be evaluated 
  + ```form-mode-2.obj```: traces and all the requested relacements will be
                           done.
- The user can create another mode. This will be achieved by modifying
  the method ```get_instructions``` in the class ```FORM```.

## What to do afer this?
- The ```form-mode-X.obj``` are the final product. After this, the files
  could be passed to other repositories for further processing:
  + latex: it takes ```form-mode-1.obj``` and generate a pdf file 
    with the squared amplitudes 
  + ParFrac: it takes ```form-mode-2.obj``` and perform partial fractions.
             This is only useful for 2-3 body phase space integration. 
  






