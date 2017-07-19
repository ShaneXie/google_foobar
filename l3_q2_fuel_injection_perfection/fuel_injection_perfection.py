import math


def answer_attempt_1(n):
    num = int(n)
    rt = math.log(num, 2)
    int_part = int(rt // 1)

    diff_1 = num - 2**int_part
    diff_2 = 2**(int_part + 1) - num

    if 2**int_part == num:
        return int(rt)
    elif diff_1 < diff_2:
        return int(int_part + diff_1)
    else:
        return int(int_part + diff_2)

cache = {1: 0, 2: 1}


def rec(n):
    if n in cache:
        return cache[n]

    if n % 2 != 0:
        cache[n] = min(rec((n + 1) / 2) + 2, rec((n - 1) / 2) + 2)
    else:
        cache[n] = rec(n / 2) + 1

    return cache[n]


def answer(n):
    # use DP
    return rec(int(n))



def test():
    # print answer(4416854613516584365132156)
    # print answer(4)
    print answer(15)
    print answer(13)

test()
