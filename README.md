# Numerical Analysis - Root Approximation Methods
Implemented PA = LU factorization, Bisection, Newton Raphson and Secant Methods.

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
  
