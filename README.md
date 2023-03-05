# Numerical Analysis - Root Approximation Methods
Implemented PA = LU factorization, Bisection, Newton Raphson, Secant Methods and Cholesky Decomposition.

In this project I will use 2 functions.Bisection, Newton Raphson and Secant Methods are implemented on the first function and the second one is used for Modified methods I created (from Bisection Newton-Raphson and Secant Methods).

### First function

![1](https://user-images.githubusercontent.com/37185221/222966961-8f43b200-0000-43fb-8aae-ed8c9784352c.PNG)

![2](https://user-images.githubusercontent.com/37185221/222967003-7aa461fa-7081-44c1-9d86-d4f934d3221f.PNG)

F(x) has 3 roots, **x0 = −1.19762, x1 = 0, x2 = 1.53013** .

```ruby
def f(x):
    result =math.e**((math.sin(x))**3) + (x**6) - (2*(x)**4) - (x**3) - 1
    return result
```
## Bisection Method

![2](https://user-images.githubusercontent.com/37185221/222967544-14ea2148-712d-47a8-a82e-0b819da32f2f.PNG)

First, for the interval `[-2,2]` we make sure that values are of opposite sign (there is at least one zero crossing within the interval).
 ```ruby
print(f(-2),' f(-2)')
print(f(2),' f(2)')
print(f(1),' f(1)')
```

![3](https://user-images.githubusercontent.com/37185221/222967845-8b3d57b5-a227-4473-b94f-4de1c213662b.PNG)
```ruby
e = 1/2 * 10**-5

def bisection_method(a,b,iterat):
    
   
    if  f(a)*f(b)<0:
        m = (a + b ) / 2
        if (np.abs(f(m)) < e):
            m= m + e
            return m, iterat

        if f(a)*f(m)<0:
            iterat=iterat+1
            return bisection_method(a,m,iterat)
        else:
            iterat=iterat+1
            return bisection_method(m,b,iterat)
```
Constant e is defined ` e = 1/2 * 10^(-k)` where k represents the desired number accuracy. In this case, `k=5`. I run bisection for intervals `[-2,1]` and `[1,2]`.
```ruby
bisection_method(-2,1,0)
```
![4](https://user-images.githubusercontent.com/37185221/222968605-a3fec52c-7d9f-49b1-9232-5a19cc2cc2ed.PNG)


```ruby
bisection_method(1,2,0)
```
![5](https://user-images.githubusercontent.com/37185221/222968716-ae213cbc-bcc6-4854-8289-583998939fba.PNG)

## Newton-Raphson Method
Newton–Raphson method starts with a single-variable function f defined for a real variable x, the function's derivative f′, and an initial guess x0 for a root of f. If the function satisfies sufficient assumptions and the initial guess is close, then:

![6](https://user-images.githubusercontent.com/37185221/222969270-9bee395a-50ac-4bfe-8cff-94e0be13fa8c.PNG)

is a better approximation of the root than x0. Geometrically, (x1, 0) is the intersection of the x-axis and the tangent of the graph of f at (x0, f(x0)): that is, the improved guess is the unique root of the linear approximation at the initial point. The process is repeated as:

![7](https://user-images.githubusercontent.com/37185221/222969318-80b3e889-fcfb-4f53-93df-eb3141e2a7db.PNG)

until a sufficiently precise value is reached. The number of correct digits roughly doubles with each step. <sub>Source: [Wiki](https://en.wikipedia.org/wiki/Newton%27s_method)</sub>

Derivative f′:
```ruby
def derf(x):
    result = 3*math.e**((math.sin(x))**3)*math.cos(x)*(math.sin(x))**2 + 6*x**5 -8*x**3 - 3*x**2
    return result
```
Method:
```ruby
def newton_raphson_method(x):
    iterat=0
    while (np.abs(f(x)) > e):
        iterat=iterat+1
        x = x - f(x)/derf(x)
        #print(x, iterat)
    return x, iterat
```
Run method for variables `x0=2, x1=-2, x2=0`.
```ruby
newton_raphson_method(2)
```
![8](https://user-images.githubusercontent.com/37185221/222969725-117e2e87-1b59-4d55-b19d-0309136517e4.PNG)

```ruby
newton_raphson_method(-2)
```
![9](https://user-images.githubusercontent.com/37185221/222969751-989f5434-eca2-4aac-889b-da76b2cde418.PNG)

```ruby
newton_raphson_method(0)
```
![10](https://user-images.githubusercontent.com/37185221/222969762-78d0c87e-7985-498d-abd1-f6c14d884f3b.PNG)

## Secant Method
For finding a zero of a function f, the secant method is defined by the recurrence relation.
![11](https://user-images.githubusercontent.com/37185221/222969838-e7f59534-f43a-4f3d-b972-4bcee785aad3.PNG)

As can be seen from this formula, two initial values x0 and x1 are required. Ideally, they should be chosen close to the desired zero.<sub>Source: [Wiki](https://en.wikipedia.org/wiki/Secant_method)</sub>

```ruby
def secant_method(a,b):
    iterat=0
    while (np.abs(f(a)- f(b)) > e):
        iterat=iterat+1
        temp=a
        a = a - (a-b) * f(a)/(f(a)-f(b))
        b= temp
    return a,iterat
  ```
  
For intervals `[-2,2]` and `[-2,0]`:
```ruby
secant_method(-2,2)
```
![12](https://user-images.githubusercontent.com/37185221/222974106-4079e373-1b1b-445b-8551-c3b6ccb021d9.PNG)

```ruby
secant_method(-2,0) 
```
![13](https://user-images.githubusercontent.com/37185221/222974927-7a682035-5165-43a3-a270-8268da7873e1.PNG)

### Second function
![14](https://user-images.githubusercontent.com/37185221/222976288-d00a1729-03fd-4c66-ae9d-900a63d569c7.PNG)

![9](https://user-images.githubusercontent.com/37185221/222976295-efb71b28-392d-4e62-ab22-b42385723276.PNG)

## Modified Newton-Raphson 
It's Newton-Raphson where the process is repeated as:

![15](https://user-images.githubusercontent.com/37185221/222976553-c8317f64-b7d6-4696-b5aa-c9ea4d5f6950.PNG)

until a sufficiently precise value is reached.

```ruby
def modified_NR(f2,derf2,dderf2,x):
    iterat=0
    while (np.abs(f2(x)) > e):
        iterat=iterat+1
        x= x - 1/((derf2(x)/f2(x)) - 1/2*(dderf2(x)/derf2(x))) #modified
    return x,iterat
```
I compare results:
```ruby
modified_NR(f2,derf2,dderf2,1)
```
![16](https://user-images.githubusercontent.com/37185221/222976729-6f275a4b-11b0-40a4-8d5a-67fda07c23f0.PNG)

```ruby
modified_NR(f2,derf2,dderf2,3)
```
![17](https://user-images.githubusercontent.com/37185221/222976758-5e410f72-bf7e-4e63-8774-e05a46147fb2.PNG)

```ruby
modified_NR(f2,derf2,dderf2,2)
```
![18](https://user-images.githubusercontent.com/37185221/222976800-d274bedc-9203-4b9d-a5b5-ccd1953a96c7.PNG)

## Modified Bisection Method
A Bisection Method where the estimation of the root is not the mean value of the interval, but a random value (each iteration).

```ruby
def modified_BM(a,b,iterat):    
   
    if  f2(a)*f2(b)<0:
        m = random.uniform(a,b) #modified
        if (np.abs(f2(m)) < e):
            m= m + e
            return m, iterat

        if f2(a)*f2(m)<0:
            iterat=iterat+1
            return modified_BM(a,m,iterat)
        else:
            iterat=iterat+1
            return modified_BM(m,b,iterat)
 ```

Testing: 
```ruby
modified_BM(0.9,2,0)
```
![20](https://user-images.githubusercontent.com/37185221/222977278-8ac15aef-05ee-4436-9f82-beb267dfaf1b.PNG)

## Modified Secant Method
A Secant method where we start from three initial values. This method is defined by the recurrence relation:

![21](https://user-images.githubusercontent.com/37185221/222977449-66c589f5-d1c3-45b0-8169-f04b47b0ec22.PNG)

where,

![22](https://user-images.githubusercontent.com/37185221/222977494-0b59835e-8132-4de0-88ac-92d315c34487.PNG)

```ruby
def modified_SM(x1,x2,x3): 
    iterat=0
    q=f2(x1)/f2(x2)
    r=f2(x3)/f2(x2)
    s=f2(x3)/f2(x1)
    while (np.abs(f2(x3)) > e):
        iterat=iterat+1
        x4 = x3 - (r*(r-q)*(x3-x2)+(1-r)*s*(x3-x1))/((q-1)*(r-1)*(s-1))
        x1=x2
        x2=x3
        x3=x4
    return x3,iterat
```
Testing:
```ruby
modified_SM(0,1,3)
```
![23](https://user-images.githubusercontent.com/37185221/222977570-f04607b2-86fb-4e81-94b8-03fd452fbb84.PNG)

## PA=LU
An LU factorization refers to the factorization of A, with proper row and/or column orderings or permutations, into two factors – a lower triangular matrix L and an upper triangular matrix U: **A=LU**  

![24](https://user-images.githubusercontent.com/37185221/222977937-1d426f1c-a2c0-451b-83ea-1935f062f4fc.PNG)

LU factorization with partial pivoting (LUP) refers often to LU factorization with row permutations only:
**PA=LU**
,where L and U are again lower and upper triangular matrices, and P is a permutation matrix, which, when left-multiplied to A, reorders the rows of A.<sub>Source: [Wiki](https://en.wikipedia.org/wiki/LU_decomposition)</sub>

First I find L, P and U matrixes.
```ruby
def find_LPU(A):
    
    L = np.zeros((A.shape[0],A.shape[1]))  #initialize with 0
    P= np.diag(np.ones(A.shape[0]))        #initialize: 1 diagonial 0 elsewhere.
    
    for j in range(A.shape[0]-1): 
            if (A[j][1]==0 or A[j][1]<A[j+1][1]):
                Pnew= np.diag(np.ones(A.shape[0]))
                A[[j+1,j],:]=A[[j,j+1],:] 
                
                Pnew[[j+1,j],:]=Pnew[[j,j+1],:] #for every line shift of A , do the same to P.
                P = np.dot(Pnew,P)              #multiple older with present P

    for i in range(0, A.shape[0]): 
       
            
        L[i][i]=1

        for j in range(i + 1, A.shape[0]): #Gauss
            if (A[i][j]): 
                l=A[j][i]/A[i][i]
                A[j] = A[j] - A[i] * (A[j][i]/A[i][i])
                L[j][i] = l
        U=A
                
    return L,P,U
    ```
   Implement PA=LU:
   ```ruby
   def PALU(A,b):
    L,P,U= find_LPU(A)
    b = np.dot(P,b) 
    print("L\n",L)
    print("P\n",P)
    print("U\n",U)
    print(b)
    y = solve(L,b) #Ly=b
    print("y= ",y) 
    x=solve(U,y)   #Ux=y
    
    return x
 ```
Test method for matrix `A=[2,1,5],[4,4,-4],[1,3,1]` and `b=[5,0,6]`.
```ruby
M= np.array([[2,1,5],[4,4,-4],[1,3,1]])
b=np.array([5,0,6])
PALU(M,b)
```
![25](https://user-images.githubusercontent.com/37185221/222978343-0eb15d67-f5f2-4644-963f-4ec329ee2a05.PNG)

## Cholesky Decomposition
s a decomposition of the form **A=LL*** where L is a lower triangular matrix with real and positive diagonal entries, and L* denotes the conjugate transpose of L.<sub>Source: [Wiki](https://en.wikipedia.org/wiki/Cholesky_decomposition)</sub>

```ruby
def cholesky(A):
    L = np.array([[0.0] * A.shape[0] for i in range(A.shape[0])])
   
    for i in range(A.shape[0]):
        for j in range(i+1):
            sum_ch = sum(L[i][k] * L[j][k] for k in range(j))
            
            if (i == j): 
                L[i][j] = math.sqrt(A[i][i] - sum_ch)
            else:
                L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum_ch))
    return L
 ```
 Test for matrix `A=[6, 3, 4, 8], [3, 6, 5, 1], [4, 5, 10, 7], [8, 1, 7, 25]` :
 ```ruby
 A = np.array([[6, 3, 4, 8], [3, 6, 5, 1], [4, 5, 10, 7], [8, 1, 7, 25]])
L = cholesky(A)
print("A:\n",A)


print("L:\n",L)
```
![26](https://user-images.githubusercontent.com/37185221/222978570-61f782bc-3906-444e-81b9-8abd32c5d489.PNG)


