def answer(xs):
    pos = [x for x in xs if x > 0]
    neg = [x for x in xs if x < 0]

    pos_count = len(pos)
    neg_count = len(neg)

    if (pos_count == 0 and neg_count == 0) or (pos_count == 0 and neg_count == 1):
        return str(xs[0])
    elif pos_count == 0:
        if neg_count % 2 == 1:
            neg.remove(max(neg))
        return str(reduce(lambda a, b: a * b, neg))
    elif neg_count == 0 or neg_count == 1:
        return str(reduce(lambda a, b: a * b, pos))
    else:
        if neg_count % 2 == 1:
            neg.remove(max(neg))
        return str(reduce(lambda a, b: a * b, pos) * reduce(lambda a, b: a * b, neg))


print answer([-1])
print answer([2])
print answer([0])
print answer([-10, -2])
print answer([-10, -2, -3])
print answer([-10, 1, 2])
print answer([-2, -3, 4, -5])
