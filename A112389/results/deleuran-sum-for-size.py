import sys
import math

S = {} # Symmetric refinements: string -> cnt
N = {} # Non-symmetric

# Set all A(X,X):
N['11'] = 22
S['11'] = 2
XX = '1'
XXs = 2
XXn = 22
for i in range(3,12):
    # An(X,X) = 46*An(X-1,X-1) + 6*2^X - 2^(X-1)
    # As(X,X) = 2 * As(X-1,X-1) = 2^(X-1)
    XX = XX + '1'
    S[XX] = XXs
    XXs = XXs * 2
    N[XX] = XXn
    XXn = 46*XXn + 6*2**i - 2**(i-1)

def setData(data):
    for i in range(0, len(data), 3):
        N[data[i]] = data[i+1] - data[i+2]
        S[data[i]] = data[i+2]

# Size 1:
setData(['1', 1, 1])

# Size 2:
# Done above

# Size 3:
setData(['21', 250, 20])

# Size 4:
setData(['31', 648, 8, '22', 10411, 49, '121', 37081, 32])

# Size 5:
setData(['41', 550, 28, '32', 148794, (443), '131', 433685, (24), '221', 1297413, (787)])

# Size 6:
setData(['51', 138, (4), '42', 849937, (473), '33', 6246077, (432)])
setData(['141', 2101339, (72), '321', 17111962, (671), '231', 41019966, (1179), '222', 43183164, (3305)])
setData(['1221', 157116243, (663)])

# Size 7:
setData(['61', 10, (4), '52', 2239070, (1788), '43', 106461697, (10551)])
setData(['421', 94955406, (6066), '241', 561350899, (15089), '331', 1358812234, (1104), '322', 561114147, (17838), '232', 3021093957, (46219), '151', 4940606, 12])
setData(['2221', 5227003593, (33392), '1321', 4581373745, (1471)])

# Size 8:
setData([
    '62', 2920534, (830),
    '53', 884147903, (5832),
    '44', 4297589646, (34099),
    '521', 245279996, (2456),
    '431', 20790340822, (23753),
    '422', 3125595194, (26862),
    '341', 41795025389, (17430),
    '332', 90630537410, (52944),
    '323', 7320657167, (14953),
    '251', 3894847047, (9174),
    '242', 84806603578, (143406),
    '161', 6059764, (12),
    '3221', 68698089712, (14219),
    '2321', 334184934526, (47632),
    '2231', 150136605052, (48678),
    '2222', 174623815718, (191947),
    '1421', 60442092848, (8871),
    '1331', 287171692047, (2640),
    '12221', 625676928843, (19191),
])

# Size 9:
setData([
    '72', 1989219, 1895,
    '63', 3968352541, 58092,
    '54', 82138898127,	281500,

    '621', 315713257,	10343,
    '531', 163360079558, (12990),
    '522', 8147612224,	74040,
    '441', 1358796413148, (525989),
    '432', 1324027972321, (901602),
    '423', 41469827815,	143968,
    '351', 647955015327, (16302),
    '342', 4999009855234, (1460677),
    '333', 2609661915535, (52782),
    '261', 15217455035,	68536,
    '252', 1221237869323, (895646),
    '171', 4014751, 0,

    '4221', 392742794892,	301318,
    '3321', 10036269263050, (59722),
    '3231', 1987600812703,	33113,
    '3222', 2312168563229,	759665,
    '2421', 8997607757089, (931275),
    '2331', 18957705069902, (119960),
    '2322', 10986279694674, (1941786),
    '2241', 1976231834547, (659723),
    '1521', 412118298729, (6758),
    '1431', 7941161106368, (37388),

    '22221', 20883741916735, (1455759),
    '13221', 17976842184698, (33957),
    '12321', 36790675675026, (39137),
])

# Size 10:
setData([
    '82', 709854, (316),
    '73', 10301630152, (21402),
    '64', 859832994275, (499397),
    '55', 3205349758318, (286406),

    '721', 212267872, (2325),
    '631', 709239437077, (122742),
    '622', 10610010722, (42938),
    '532', 10198551751032, 592088,
    '523', 110432745036, (58784),
    '541', 23168524352411, 435708,
    '451', 41531542406815, 772386,
    '442', 144735111618598, (5784742),
    '433', 37566339738080, (1069641),
    '424', 241236702180, (221465),
    '361', 5711086649169, (112022),
    '352', 134764333145996, 1639902,
    '343', 262440584015903, (1688509),
    '271', 35758538164, (24913),
    '262', 10134629875966, (1466770),
    '181', 1421072, 0,

    '5221', 1064278709384, (55376),
    '4322', 4914171473466769, (38340865),
    '4321', 147793134818751, (808943),
    '4231', 11653960252958, (414800),
    '4222', 13378142987817, (1629981),
    '3421', 528069494287014, 729975,
    '3331', 547495815712759, (123794),
    '3322', 331549223161406, (2210342),
    '3241', 26468746650129, (206075),
    '3232', 147000420605317, (1060478),
    '3223', 30853217686804, (303826),
    '2521', 126768194057206, 779906,
    '2431', 936478355031379, 3294187,
    '2422', 294418057243489, (7325886),
    '2341', 500091779357026, (1310278),
    '2332', 1250049347446753, (4902183),
    '2251', 13526583972859, (398785),
    '1621', 1592586147307, (21527),
    '1531', 116300201229509, 34810,
    '1441', 419826910043616, 493182,

    '32221', 277488918507907, (421588),
    '23221', 1305158898588543, (1139342),
    '22321', 1209535848675777, (2034360),
    '22231', 602787318883898, (2094327),
    '22222', 697608586669144, (10421527),
    '14221', 238416260244308, (309230),
    '13321', 2079934426148637, (128102),
    '13231', 518058446706002, (74915),
    '12421', 952602938632840, (359320),

    '122221', 2488886491814997, (628498),
])

# Size 11:
setData([
    '92', 129568, (552),
    '83', 16200206750, (112636),
    '74', 5357035940501, (2290271),
    '65', 66349485360974, 7324963,

    '821', 75044114, (4916),
    '731', 1799186992768, 44114,
    '722', 7211055824, (80482),
    '641', 226868353416156, 5961252,
    '632', 43942851658601, 4869223,
    '623', 147204185237, (260083),
    '551', 1338785905226577, (718180),
    '542', 2315339676520986, (32326444),
    '533', 289481658870354, (632360),
    '524', 662563743656, (629320),
    '461', 708637636378386, 11167524,
    '452', 7017458196473746, (59483570),
    '443', 7150024883019288, (44315210),
    '434', 541700127346014, (17080083),
    '371', 30937971078448, 61287,
    '362', 2077242826111952, (27498922),
    '281', 52647227697, (118808),
    '272', 52338565807622, 5846935,
    '191', 258584, 0,

    '6221', 1464493253086, (667311),
    '5321', 1155060203595226, (588844),
    '5231', 32708017336078, (132016),
    '5222', 36851077736763, (3166928),
    '4331', 7985866751161543,	2371105,
    '4322', 4914171473466769, (38340865),
    '4241', 158892437059818, (6283476),
    '4232', 879794762964609, (17193399),
    '4223', 180217829542618, (6905133),
    '3521', 13876950263616988, (1367292),
    '3422', 17332556350873758, (64238315),
    '3341', 14226359474548568, (1612872),
    '3332', 36288556463260004, (5580696),
    '3323', 4472233899139020, 1348484,
    '3242', 3983141281731531, (21660689),
    '3251', 183657614407425, (164712),
    '2621', 1039388269581582, (8373861),
    '2522', 4149870746820242, (38305026),
    '2351', 7154845370835644, (1495234),
    '2261', 52566014594439, 3074177,
    '1721', 3710232065761, 8476,
    '1631', 991756464038495, (232826),

    '42221', 1619895602468513, (13822233),
    '33221', 39335472994895589, (1402284),    
    '32231', 8071935524995532, (730718),
    '32222', 9286460454529759, (33185404), # Old Lemma 2 algorithm counted +1 by mistake!
    '24221', 34964858265262896, (44673646),
    '23321', 136668251667112320, (4613691),
    '23231', 37712858319195719, (2418614), # Old Lemma 2 algorithm counted +1 by mistake!
    '23222', 43743183773027066, (84481663), # Old Lemma 2 algorithm counted incorectly!
    '22322', 39797545797160980, (81467056),
    '22331', 68666008843350491, (5033618),
    '22421', 31179019407650704, (41370675),
    '22241', 8042576327798896, (29049583),
    '15221', 1654910007480680, (200521),
    '14321', 54989029357667553, (1080236),
    '14231', 6951175887318281, (519900),
    '13331', 112790108951168181, (284498),
    '12521', 13147170177676408, (362230),

    '222221', 83131865065198060, (64343390), # Old Lemma 2 algorithm counted incorectly!
    '132221', 71849872746311779, (1046044), # Old Lemma 2 algorithm counted incorectly!
    '123221', 143351914222644371, (1028025),
])

# Impossible constructions:
setData(['91', 0, 0, '81', 0, 0, '71', 0, 0])


def get(prefix, A):
    if prefix in A:
        return A[prefix]
    if prefix[::-1] in A:
        return A[prefix[::-1]]
    if len(prefix) < 3:
        return None
    for i in range(1, len(prefix)-1):
        if prefix[i] == '1':
            prefixL = prefix[0:i+1:]
            left = get(prefixL, A)
            prefixR = prefix[i::]
            right = get(prefixR, A)
            if left != None and right != None:
                return [prefixL, prefixR]

def getS(prefix):
    x = get(prefix, S)
    if x == None or type(x) is int:
        return x
    cs = get(x[0], S)
    ds = get(x[1], S)
    S[prefix] = cs * ds
    return cs * ds

def getN(prefix):
    x = get(prefix, N)
    if x == None or type(x) is int:
        return x
    cn = get(x[0], N)
    cs = get(x[0], S)
    dn = get(x[1], N)
    ds = get(x[1], S)
    d = dn+ds
    ret = (dn+d)*cn + d*cs - cs*ds
    N[prefix] = ret
    return ret


printed = {}
cntRefinements = 0
cntLemma1 = 0
def printXY(X, Y, prefix, rem, actuallyPrint):
    global cntRefinements, cntLemma1

    if rem == 0 and len(prefix) == Y:
        s = getS(prefix)
        n = getN(prefix)
        sum = [0,0] if None in [s,n] else [n,s]
        if prefix in printed:
            return sum # Already printed: Should still be counted
        printed[prefix] = True
        printed[prefix[::-1]] = True
        if n == 0 and s == 0 or not actuallyPrint:
            return sum
        cntRefinements = cntRefinements + 1
        if '1' in prefix[1:-1]:
            cntLemma1 = cntLemma1 + 1
        if n != None and s != None:
            print('  <' + prefix + '>', s + getN(prefix), '(' + str(s) + ')')
        else:
            print('  <' + prefix + '> #')
        return sum
    if len(prefix) >= Y:
        return (0,0) # Too tall!
    sum = [0,0]
    for c in range(min(9, rem), 0, -1):
        add = printXY(X, Y, prefix + str(c), rem-c, actuallyPrint)
        sum[0] = sum[0] + add[0]
        sum[1] = sum[1] + add[1]
    return sum

to = int(sys.argv[1])
print()
print('Refinements of size', to)
for X in range(2,to+1):
    sumN = 0
    sumS = 0
    for Y in range(2,X+1):
        if X == to:
            print(' Height', Y)
        [n,s] = printXY(X, Y, '', X, X == to)
        sumN = sumN + n
        sumS = sumS + s
        if X == to:
            print('   SUM', n+s, '(' + str(s) + ')')
    if X == to:
        print('TOTAL', sumN+sumS, '(' + str(sumS) + ')')
print()
#print('Number of refinements:', cntRefinements)
#print('Found using Lemma 1:', cntLemma1)
#print()


# writeHeatmapFile() writes the heatmap.htm file which shows colors
# of refinements based on the amount of non-symmetric counts
def writeHeatmapRow(f, size, height, prefix, rem, sizeMin, sizeMax):
    if rem == 0 and len(prefix) == height:
        n = getN(prefix)
        if prefix in printed:
            return # Already printed: Should still be counted
        printed[prefix] = True
        printed[prefix[::-1]] = True
        if n == 0:
            return # Impossible refinement
        if '1' in prefix[1:-1]:
            return # Lemma 1
        if n != None and sizeMax != sizeMin:
            nom = n-sizeMin
            denom = sizeMax-sizeMin
            #gradient = 1 - (math.log(1+nom)/math.log(1+denom))
            #gradient = 1 - (nom*nom/(denom*denom))
            gradient = 1 - nom/denom
            gradient = 100 * gradient * gradient
            f.write('<td style="background-color:hsl(' + str(gradient) +',65%,50%);">')
        else:
            f.write('<td>')
        nn = '%.0E' % n if n != None else ''
        f.write('&lt;' + prefix + '&gt;<br/>' + nn + '</td>')
        return
    if len(prefix) >= height or rem < 1:
        return
    for c in range(min(9,rem), 0, -1): # No 2+-digit rows
        writeHeatmapRow(f, size, height, prefix + str(c), rem-c, sizeMin, sizeMax)
def writeHeatmapFile():
    with open('heatmap.htm', 'w') as f:
        f.write('<html>\n')
        f.write(' <head>\n')
        f.write('  <style>\n')
        f.write('   td {text-align: center; padding: 2px;}\n')
        f.write('  </style>\n')
        f.write(' </head>\n')
        f.write(' <body>\n')
        for size in range(8,to+1):
            # Print:
            f.write('  <h3>Size ' + str(size) + '</h3>\n')
            f.write('  <table>\n')
            for height in range(2,size+1):
                # Find min and max:
                sizeMin = 1e20
                sizeMax = 0
                for key in N:
                    if len(key) == height and sum([int(x) for x in key]) == size and not '1' in key[1:-1]:
                        sizeMin = min(sizeMin, N[key])
                        sizeMax = max(sizeMax, N[key])
                f.write('   <tr>\n')
                writeHeatmapRow(f, size, height, '', size, sizeMin, sizeMax)
                f.write('   </tr>\n')
            f.write('  </table>\n')
        f.write(' </body>\n')
        f.write('</html>\n')
printed = {} # Reset
writeHeatmapFile()


# setupKnownCounts() computes and prints cross check values for "Combination::setupKnownCounts()"
def setupKnownCounts():
    global S, N

    for k in N:
        if '1' in k[1:-1]:
            continue
        n = N[k]
        if n <= 1:
            continue
        s = S[k]
        a = str(n+s)
        s90 = '122' if k == '44' else '0'
        print(' '*4 + 'm[' + k + '] = Counts(' + a + ', ' + str(s) + ', ' + s90 + ');')
    pass
setupKnownCounts()


# Code below recreates Table 7 from Eilers (2016):
def prep():
    return [[0 for x in range(to)] for y in range(to)]
f = [0] * to # c with line on top in Table 7
f180 = [0] * to
C = prep() # c(n,m) with bold c
C90 = prep()
if to >= 8:
    C90[7][3] = 244
C180 = prep()

def fillC(n, prefix, rem):
    if rem == 0 and len(prefix) > 0:
        if prefix[0] == '1':
            # f and f180:
            fprefix = prefix + '1'
            N = getN(fprefix)
            S = getS(fprefix)
            if N is not None:
                f[n-1] = f[n-1] + 2 * N + S
                f180[n-1] = f180[n-1] + getS(fprefix)
        # C and C180:
        firstLayer = int(prefix[0])
        N = getN(prefix)
        S = getS(prefix)        
        if N is not None:
            C[n-1][firstLayer-1] = C[n-1][firstLayer-1] + (2 * N + S) * firstLayer
            C180[n-1][firstLayer-1] = C180[n-1][firstLayer-1] + S * firstLayer
        return
    for m in range(1 if prefix == '' else 2, rem+1):
        fillC(n+m, prefix + str(m), rem-m)

def multf(remainingCount, remainingN, product, ret, go180):
    if remainingCount == 0:
        if remainingN != 0:
            return
        ret.append(product);
        return
    for i in range(1, remainingN+1):
        multf(remainingCount-1, remainingN-i, product * (f180[i-1] if go180 else f[i-1]), ret, go180)

def eilers():
    for n in range(1, to+1):
        fillC(0, '', n)

    print('n\tf(n)\tf180(n)\t' + '\t'.join(["C(n,{0})\tC180(n,{0})".format(x) for x in range(1,to+1)]))
    for n in range(0, to):
        printList = [n+1,f[n],f180[n]]
        for m in range(0, to):
            printList.append(C[n][m])
            printList.append(C180[n][m])
        print('\t'.join([str(x) for x in printList]))

    # Compute a(n) using 2.3 from Eilers (2016):
    for n in range(1, to+1):
        an = 0
        # First row sum:
        for m in range(2, n+1):
            an = an + (C[n-1][m-1] + C180[n-1][m-1] + 2 * C90[n-1][m-1]) / (2 * m)
        # Second row first sum:
        suml = 0
        for l in range(n): # l = 0...n-1
            for m1 in range(1,n+1): # m1 = 1..n
                for m2 in range(1,n+2-m1): # m2 = 1..
                    sums = []
                    multf(l, n+1-m1-m2, 1, sums, False)
                    for sum in sums:
                        suml = suml + C[m1-1][0] * C[m2-1][0] * sum
                    sums180 = []
                    multf(l, n+1-m1-m2, 1, sums180, True)
                    for sum in sums180:
                        suml = suml + C180[m1-1][0] * C180[m2-1][0] * sum
        an = an + suml / 2
        print('a(' + str(n) + ') =', an)
# eilers() # Run this line to compute a(n) and table 2.3 from Eilers (2016)
