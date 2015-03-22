# PyPCS
Python framework that computes several kind of parametric curves and surfaces.

## About
This framework has been developed as part of my bachelor's degree thesis work. It uses [NumPy](http://www.numpy.org/)
and [Matplotlib](http://matplotlib.org/).
The framework contains code to compute:

         | Curve  | Surface  
----------|--------|-----------
 Bézier  |  Yes   |   Yes    
 Natural Spline |   Yes | No
 Hermite Spline | Yes | No
 Cardinal Spline | Yes | No
 B-Spline | Yes | Yes
 NURBS    | Yes | Yes

It is also possible to perform curve and surface fitting using NURBS.

## Usage
To create a curve, you have to provide a list of control points (or the name of a file that contains them) and the number of points of the curve you want to compute. For example, for a Bézier curve

```python
    curve = Bezier([(0, 5), (5, 8), (10, 8), (15, 3)], 200)
    curve.calculate()
```
