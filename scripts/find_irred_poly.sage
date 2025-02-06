
import itertools

K.<x> = GF(16)

coefs = [0,1,x,x**2,x**3]

score = {}
for k in K:
    score[k] = 0;

for k in coefs:
    score[k] = 1

score[0] = 4
score[1] = 2

PR.<z> = PolynomialRing(K)

for m in [78,64,108,142]:

    tries = 0
    best_score = -1
    for tail in itertools.product(K, repeat = 4):
        P = tail[0] + tail[1]*z + tail[2]*z**2 + tail[3]*z**3 + z**m  

        s= sum([score[k] for k in tail])
        if tail[-1] == 0:
            s += 2
        if s <= best_score:
            continue

        tries += 1

        if P.is_irreducible():
            best_score = s
            print("found: ", P, " score: ", best_score)
            
    print("tries:", tries)