import sympy as sy
import numpy as np
from matplotlib import pyplot as plt
from math import *
import sys
import pygame
from pygame.locals import *

word=['','강산을 강염기로 적정','강염기를 강산으로 적정','약산을 강염기로 적정','약염기을 강산으로 적정']


def resultprint(eq_volume,x,y,flag):
    global word
    print()
    print()
    str1="<"+word[flag]+"> 계산 결과"
    print(str1)
    str2="- 당량 부피: "+str(eq_volume)+"mL"
    print(str2)
    str3="- v=0mL에서의 pH: "+str(round(y[0],2))
    print(str3)
    str4="- v="+str(round(x[len(x)//4],2))+"mL(당량점/2) 에서의 pH: "+str(round(y[len(x)//4],2))
    print(str4)
    str5="- v="+str(round(x[len(x)//2],2))+"mL(당량점) 에서의 pH: "+str(round(y[len(x)//2],2))
    print(str5)
    str6="- v="+str(round(x[len(x)//4*3],2))+"mL(당량점*3/2) 에서의 pH: "+str(round(y[len(x)//4*3],2))
    print(str6)

def calculator(ca,cb,ka,kw=10**-14):
    x=sy.Symbol('x')
    equation=x**3+(cb+ka)*x**2-(kw+ca*ka)*x-ka*kw
    #print(equation)
    ans=sy.solveset(equation,x)
    for i in ans:
        if i>0:
            return -log10(i)
            #return round(-log10(i),2)
    #print(ans)

def strongstrong(acid_volume,acid_con,acid_n,base_volume,base_con,base_n):
    con=(acid_volume*acid_con*acid_n-base_volume*base_con*base_n)/(acid_volume+base_volume)
    if con>0:       #산성0
        return -log10(con)
    elif con==0:    #중성
        return 7
    else:           #염기성
        flag=3
        return 14+log10(-con)

def strong(init_volume,init_con,init_n,add_con,add_n,flag):
    #초기 부피, 초기 농도, 초기 가수, 첨가 농도
    #flag : 1(산을 염기로 적정), 2(염기를 산으로 적정)
    #print(init_volume,init_con,init_n,add_con,add_n,flag)
    eq_volume=round(init_volume*init_con*init_n/(add_con*add_n),2)
    #print(eq_volume)
    x=np.arange(0,eq_volume*2,0.1)
    y=[]
    for i in range(len(x)):
        if flag==1:
            y.append(strongstrong(init_volume,init_con,init_n,i*0.1,add_con,add_n))
        elif flag==2:
            y.append(strongstrong(i*0.1,add_con,add_n,init_volume,init_con,init_n))
    y=np.array(y)
    return x,y,eq_volume

def weakacid(init_volume,init_con,ka,add_con,add_n,flag):
    eq_volume=eq_volume=round(init_volume*init_con/(add_con*add_n),2)
    #print(eq_volume)
    x=np.arange(0,eq_volume*2,0.1)
    y=[]
    print()
    print("진행 상황 [ ",end="")
    for i in range(len(x)):
        ca=float((init_volume*init_con-add_con*0.1*i)/(init_volume+0.1*i))
        if ca>0:
            cb=float((add_con*0.1*i)/(init_volume+0.1*i))
            y.append(calculator(ca,cb,ka))
            if (i+1)%(len(x)//20)==0:
                print('#',end="")
        else:
            y.append(14+log10(-ca))
        if y[-1]<y[0]:
            y[-1]=calculator(0,cb,ka)
    print(" ]")
    y=np.array(y)
    #print(x,y)
    return x,y,eq_volume

def weakbase(init_volume,init_con,kb,add_con,add_n,flag):
    ka=1e-14/kb
    eq_volume=eq_volume=round(init_volume*init_con/(add_con*add_n),2)
    #print(eq_volume)
    x=np.arange(0,eq_volume*2,0.1)
    y=[]
    print()
    print("진행 상황 [ ",end="")
    for i in range(len(x)):
        cb=float((init_volume*init_con-add_con*0.1*i)/(init_volume+0.1*i))
        if cb>0:
            ca=float((add_con*0.1*i)/(init_volume+0.1*i))
            y.append(calculator(ca,cb,ka))
            if (i+1)%(len(x)//20)==0:
                print('#',end="")
        else:
            y.append(-log10(-cb))
        if y[-1]>y[0]:
            y[-1]=calculator(ca,0,ka)
    print(" ]")
    y=np.array(y)
    #print(x,y)
    return x,y,eq_volume

def strongstart(flag):
    flag=int(flag)
    keyword=['','강산','강염기']
    str1='적정하고자 하는 '+keyword[flag]+'의 부피 입력(단위 mL): '
    str2='적정하고자 하는 '+keyword[flag]+'의 농도 입력(단위 M): '
    str3='적정하고자 하는 '+keyword[flag]+'의 가수 입력: '
    str4='첨가하는 '+keyword[3-flag]+'의 농도 입력(단위 M): '
    str5='첨가하는 '+keyword[3-flag]+'의 가수 입력 : '
    print()
    init_volume=float(input(str1)) 
    init_con=float(input(str2))
    init_n=int(input(str3))
    add_con=float(input(str4))
    add_n=int(input(str5))
    ret=strong(init_volume,init_con,init_n,add_con,add_n,flag)
    x=ret[0]
    y=ret[1]
    eq_volume=ret[2]
    #print(x,y)
    resultprint(eq_volume,x,y,flag)
    plt.title("titration curve")
    plt.plot(x,y)
    plt.show()

def weakacidstart(flag):
    flag=int(flag)
    str1='적정하고자 하는 약산의 ka 입력(예시: 8x10^-7이면 8e-7): '
    str2='적정하고자 하는 약산의 부피 입력(단위 mL): '
    str3='적정하고자 하는 약산의 농도 입력(단위 M): '
    str4='첨가하는 강염기의 농도 입력(단위 M): '
    str5='첨가하는 강염기의 가수 입력: '
    print()
    ka=float(input(str1))
    init_volume=float(input(str2))
    init_con=float(input(str3))
    add_con=float(input(str4))
    add_n=int(input(str5))
    ret=weakacid(init_volume,init_con,ka,add_con,add_n,flag)
    x=ret[0]
    y=ret[1]
    eq_volume=ret[2]
    #print(x,y)
    resultprint(eq_volume,x,y,flag)
    plt.title("titration curve")
    plt.plot(x,y)
    plt.show()

def weakbasestart(flag):
    flag=int(flag)
    str1='적정하고자 하는 약역기의 kb 입력(예시: 8x10^-7이면 8e-7): '
    str2='적정하고자 하는 약염기의 부피 입력(단위 mL): '
    str3='적정하고자 하는 약염기의 농도 입력(단위 M): '
    str4='첨가하는 강산의 농도 입력(단위 M): '
    str5='첨가하는 강산의 가수 입력: '
    print()
    kb=float(input(str1))
    init_volume=float(input(str2))
    init_con=float(input(str3))
    add_con=float(input(str4))
    add_n=int(input(str5))
    ret=weakbase(init_volume,init_con,kb,add_con,add_n,flag)
    x=ret[0]
    y=ret[1]
    eq_volume=ret[2]
    #print(x,y)
    resultprint(eq_volume,x,y,flag)
    plt.title("titration curve")
    plt.plot(x,y)
    plt.show()

def start():
    print('<적정 곡선>')
    print('1. 강산을 강염기로 적정')
    print('2. 강염기를 강산으로 적정')
    print('3. 약산을 강염기로 적정')
    print('4. 약염기을 강산으로 적정')
    flag=input('*적정 유형 선택: ')
    if flag=='1' or flag=='2':
        strongstart(flag)
    elif flag=='3':
        weakacidstart(flag)
    elif flag=='4':
        weakbasestart(flag)
    

start()
