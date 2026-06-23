result = {'low': 25, 'normal': 25, 'high': 25}
ecuts = [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]

def generate_insert_ecuts(ecuts_old, hints):
    ecuts_new = []
    ecuts_old = sorted(ecuts_old)
    for label, hint in hints.items():
        idx = ecuts_old.index(hint)
        if idx > 0:
            tmp = [i for i in range(ecuts_old[idx - 1], ecuts_old[idx])]
            tmp.pop(0)
            ecuts_new += tmp
    ecuts_new = sorted(set(ecuts_new))
    if len(ecuts_new) == 0:
        ecuts_new = [max(5, ecuts_old[0]-5)]
    print(ecuts_new)
    return ecuts_new

result = generate_insert_ecuts(ecuts, result)
print(result)
