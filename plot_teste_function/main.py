# -*- coding: utf-8 -*-
#===============================================================================
import matplotlib.pyplot as plt
from function_to_be_bounded import *
from parametros import *
#===============================================================================
def main():
    print("mu= ",mu)
    print("r = ", r)
    print("I = ", I)
    print("mu*alpha = ", mu*alpha())
    list_image = []
    for theta in domain_theta:
        image = function_to_be_bounded(theta)
        list_image.append(image)
    critical_1 = 0
    critical_2 = np.pi*(1-r*I)
    list_image = np.array(list_image)
    plt.figure()
    plt.plot(domain_theta ,list_image)
    plt.plot(critical_1 , function_to_be_bounded(critical_1), '.', color = "red")
    plt.plot(critical_2 , function_to_be_bounded(critical_2), '.', color = 'green')
    plt.show()


if __name__ == '__main__':
    main()
