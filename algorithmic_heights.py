def fibonacci_numbers(n):
    if n == 0:
        return 0
    elif n == 1:
        return 1
    return fibonacci_numbers(n-1) + fibonacci_numbers(n-2)


if __name__ == '__main__':
    print(fibonacci_numbers(23))