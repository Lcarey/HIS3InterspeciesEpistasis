from to_import import * 

def func_find_min_s12(data,data_sigma,k):
    #print k
    s12_left=0
    s12_right=1.28
    f_s12_left=func_fit_single_exponents(data,data_sigma,s12_left,k,False)
    f_s12_right=func_fit_single_exponents(data,data_sigma,s12_right,k,False)
    s12_center=(s12_left+s12_right)*0.5
    f_s12_center=func_fit_single_exponents(data,data_sigma,s12_center,k,False)
        
    for i in range(7):
        #print s12_left,s12_center,s12_right,f_s12_left,f_s12_center,f_s12_right
        if ((f_s12_center<f_s12_left)&(f_s12_center<f_s12_right)):
            s12_left=0.5*(s12_left+s12_center)
            s12_right=0.5*(s12_right+s12_center)
            f_s12_left=func_fit_single_exponents(data,data_sigma,s12_left,k,False)
            f_s12_right=func_fit_single_exponents(data,data_sigma,s12_right,k,False)
        elif ((f_s12_center>=f_s12_left)&(f_s12_center<f_s12_right)):
            s12_right=s12_center
            f_s12_right=f_s12_center
            s12_center=(s12_left+s12_right)*0.5
            f_s12_center=func_fit_single_exponents(data,data_sigma,s12_center,k,False)
        elif ((f_s12_center>=f_s12_right)&(f_s12_center<f_s12_left)):
            s12_left=s12_center
            f_s12_left=f_s12_center
            s12_center=(s12_left+s12_right)*0.5
            f_s12_center=func_fit_single_exponents(data,data_sigma,s12_center,k,False)
        else:
            #print "stop"
            break
    #print f_s12_center,s12_center
    return f_s12_center,s12_center

def func_find_min_k(data,data_sigma):
    k_left=1.2
    k_right=14
    [f_k_left,s12_k_left]=func_find_min_s12(data,data_sigma,k_left)
    [f_k_right,s12_k_right]=func_find_min_s12(data,data_sigma,k_right)
    k_center=(k_left+k_right)*0.5
    [f_k_center,s12_k_center]=func_find_min_s12(data,data_sigma,k_center)
    for i in range(7):
        print k_left,k_center,k_right,f_k_left,f_k_center,f_k_right
        if ((f_k_center<f_k_left)&(f_k_center<f_k_right)):
            k_left=0.5*(k_left+k_center)
            k_right=0.5*(k_right+k_center)
            [f_k_left,s12_k_left]=func_find_min_s12(data,data_sigma,k_left)
            [f_k_right,s12_k_right]=func_find_min_s12(data,data_sigma,k_right)
        elif ((f_k_center>=f_k_left)&(f_k_center<f_k_right)):
            k_right=k_center
            f_k_right=f_k_center
            k_center=(k_left+k_right)*0.5
            [f_k_center,s12_k_center]=func_find_min_s12(data,data_sigma,k_center)
        elif ((f_k_center>=f_k_right)&(f_k_center<f_k_left)):
            k_left=k_center
            f_k_left=f_k_center
            k_center=(k_left+k_right)*0.5
            [f_k_center,s12_k_center]=func_find_min_s12(data,data_sigma,k_center)
        else:
            #print k_left,k_center,k_right,f_k_left,f_k_center,f_k_right
            #print "stop"
            break
    print f_k_center,k_center,s12_k_center
    return f_k_center,k_center,s12_k_center

