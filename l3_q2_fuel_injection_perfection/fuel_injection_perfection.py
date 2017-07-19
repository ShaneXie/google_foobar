import math


def answer(n):
    num = int(n)
    rt = math.sqrt(num)
    int_part = int(rt // 1)

    print int_part

    diff_1 = num - 2**int_part
    diff_2 = 2**(int_part + 1) - num

    if 2**rt == num:
        return int(rt)
    elif diff_1 < diff_2:
        return int(int_part + diff_1)
    else:
        return int(int_part + diff_2)


def test():
    print answer(4416854613516584365132156)
    print answer(15)

test()
