import sys




def findseq(seqsize, searchlength):

    #print("Looking for",seqsize, " in ", searchlength, "digits of Pi")
    digitindex=0
    collision = dict()
    srcfile=open("pi-billion.txt","rb")
    # We will read the first two bytes(which have "3.1") and skip those digits.
    skip=srcfile.read(2)
    current=bytearray(srcfile.read(seqsize))
    collision.update({bytes(current):digitindex})

    while (digitindex < searchlength):
        digitindex = digitindex+1
        for x in range(seqsize-1):
            current[x]=current[x+1]
        current[(seqsize-1)] = bytearray(srcfile.read(1))[0]
        r=collision.get(bytes(current),-1)
        if (r == -1):
            collision.update({bytes(current):digitindex})
        else:
            print("Collision found!",current,digitindex+1,r+1)
            break

    print(len(collision),"digits searched")
    srcfile.close()
    collision.clear()


maxdigits = int(sys.argv[1])
print("Seaching using",maxdigits, "digits of Pi")
for z in range(3,20):
    findseq(z,maxdigits)

