When I run `conda update conda`, I get the following error:

```
Solving environment: done

## Package Plan ##

  environment location: /home/ivan16/anaconda3

  added / updated specs: 
    - conda


The following packages will be UPDATED:

    conda: 4.5.4-py36_0 --> 4.5.11-py36_0

Proceed ([y]/n)? y

Preparing transaction: done
Verifying transaction: done
Executing transaction: failed
ERROR conda.core.link:_execute(502): An error occurred while uninstalling package 'defaults::conda-4.5.4-py36_0'.
PermissionError(13, 'Permission denied')
Attempting to roll back.

Rolling back transaction: done

PermissionError(13, 'Permission denied')
```

Then just run:

```
sudo chown -R ivan16 anaconda3 
```

here ivan16 is the user name. And then you can update successfully:

```
Solving environment: done

## Package Plan ##

  environment location: /home/ivan16/anaconda3

  added / updated specs: 
    - conda


The following packages will be UPDATED:

    conda: 4.5.4-py36_0 --> 4.5.11-py36_0

Proceed ([y]/n)? y

Preparing transaction: done
Verifying transaction: done
Executing transaction: done
```
