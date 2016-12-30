import numpy as np
import h5py as h5
import sys

def compareLists( l1, l2 ):
    if ( len(l1) != len(l2) ):
        return False

    for val in l1:
        if ( not val in l2 ):
            return False
    return True

def main( argv ):
    file1 = ""
    file2 = ""
    for arg in argv:
        if ( arg.find("--help") != -1 ):
            print ("Usage: python compareH5Files.py --file1=<h5file> --file2=<h5file>")
        elif ( arg.find("--file1=") != -1 ):
            file1 = arg.split("--file1=")[1]
        elif ( arg.find("--file2=") != -1 ):
            file2 = arg.split("--file2=")[1]

    # Check that both filenames are given
    if ( file1 == "" ):
        print ("File 1 is not given!")
        return
    elif ( file2 == "" ):
        print ("File 2 is not given!")
        return

    # Compare that the same datasets are present
    with h5.File(file1, 'r') as hf1:
        dsets1 = list(hf1.keys())
        with h5.File(file2, 'r') as hf2:
            dsets2 = list( hf2.keys() )


            keysAreEqual = compareLists( dsets1, dsets2 )
            if ( keysAreEqual ):
                print ("\033[0;32mKey test passed\033[0m")
            else:
                print ("\033[0;31mThe files has not the same keys\033\[0m")
                return

            # Compare datasets
            for key in dsets1:
                val1 = np.array( hf1.get(key) )
                val2 = np.array( hf2.get(key) )
                equalSize = True
                if ( len(val1.shape) != len(val2.shape) ):
                    equalSize = False
                else:
                    for i in range(0, len(val1.shape) ):
                        if ( val1.shape[i] != val2.shape[i] ):
                            equalSize = False
                            break

                if ( not equalSize ):
                     print ( "\033[0;31mArrays %s have different dimensions\033[0m"%(str(key)) )
                     continue

                if ( np.allclose(val1, val2, rtol=1E-3) ):
                    print ( "\033[0;32mArrays %s passed\033[0m"%(str(key)) )
                else:
                    print ( "\033[0;31mArrays %s failed\033[0m"%(str(key)) )

if __name__ == "__main__":
    main( sys.argv[1:] )
