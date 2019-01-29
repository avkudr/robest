Fitting your model
----------------------------------

Goal
~~~~

   In this tutorial we will you explain how to add a new fitting model, using the example of creating `circleFitting`.

Task details
~~~~~~~~~~~~

   The input of the algorithm is a set of 2D points, which are the measured coordinates of a circle.
   This set can be very noisy and distorted (outliers are present). So, we need to implement the main body
   of the `circleFitting` class.

Code
~~~~
    
    The tutorial code's is shown lines below.

    - Header file code:    

    .. code-block:: c++
        
       /**
       *  @brief Circle fitting to a set of 2D points
       *
       *  @author  Andrey Kudryavtsev (avkudr.github.io)
       *  @author  Rahima Djahel (github:rahma24000)
       *  @date    20/03/2018
       *  @version 1.0
       */

       #include <list>
       #include <vector>
       #include <cmath>
       #include <iostream>

       #include "robust_estim.hpp"

       struct Point2d{
           double x;
           double y;
       };

       typedef std::vector<Point2d> Point2Dvector;

       class CircleFittingProblem : public robest::EstimationProblem{

       public:
           CircleFittingProblem();
           ~CircleFittingProblem();

           void setData(std::vector<double> & x, std::vector<double> & y);

           double estimErrorForSample(int i);
           void   estimModelFromSamples(const std::vector<int> & samplesIdx);

           int getTotalNbSamples() const{
	       return (int) points.size();
           }

           void getResult(double & res_cx, double & res_cy, double & res_r) const{
	       res_cx = this->cx;
	       res_cy = this->cy;
	       res_r  = this->r;
           }

           bool isDegenerate(const std::vector<int> & samplesIdx);

       private:

           Point2Dvector points; // input data
           double cx;
           double cy;
           double r; //radius
       };


    - Source file code:

    .. code-block:: c++

       #include "CircleFitting.hpp"

       CircleFittingProblem::CircleFittingProblem(){
           setNbParams(3);
           setNbMinSamples(3);
       }

       CircleFittingProblem::~CircleFittingProblem(){

       }

       void CircleFittingProblem::setData(std::vector<double> & x, std::vector<double> & y)
       {
           points.clear();
           for (int i = 0; i < x.size(); i++){
               Point2d data;
               data.x=x[i];
               data.y=y[i];
               points.push_back(data);
           }
       }

       inline double CircleFittingProblem::estimErrorForSample(int i)
       {
           // distance circle-point = abs(<distance point-center> - radius)
           Point2d & p = points[i];
           return std::abs(sqrt((p.x-cx)*(p.x-cx)+(p.y-cy)*(p.y-cy)) - r);
       }

       inline void CircleFittingProblem::estimModelFromSamples(const std::vector<int> & samplesIdx){
           if( !isDegenerate(samplesIdx)){
               Point2d & P = points[samplesIdx[0]];
               Point2d & V = points[samplesIdx[1]];
               Point2d & K = points[samplesIdx[2]];
               //calculation of the coefficients of the mediating straight lines
               double a = -(V.x - P.x)/(V.y - P.y);
               double b = (V.x * V.x - P.x * P.x + V.y * V.y - P.y * P.y)/(2* (V.y - P.y));
               double c = -(K.x - V.x)/(K.y - V.y);
               double d = (K.x * K.x - V.x * V.x + K.y * K.y - V.y * V.y)/(2* (K.y - V.y));

               //calculate the coordinates of the center of the circle O(A,B)
               cx = (b-d)/(c-a);
               cy = a*cx + b;

               //calculate the radius of a circle
               r = sqrt((P.x - cx)*(P.x - cx)+(P.y - cy)*(P.y - cy));
           }
       }

       inline bool CircleFittingProblem::isDegenerate(const std::vector<int> & samplesIdx)
       {
           Point2d & P = points[samplesIdx[0]];
           Point2d & V = points[samplesIdx[1]];
           Point2d & K = points[samplesIdx[2]];

           // verify that points P, V and K are not at the line -> verify that PV and PK are colinear:

           //1. calculate the directing coefficient of the line PV
           double f = (V.y-P.y)/(V.x-P.x);

           //2. calculate the directing coefficient of the line PK
           double h = (K.y-P.y)/(K.x-P.x);

           //3. PV and PK Are colineaire if and only if f = h
           return ( f - h < 1e-3 );
       }

    - Test file code:

    .. code-block:: c++

       #include "gtest/gtest.h"

       #include <random>

       #include "CircleFitting/CircleFitting.hpp"

       //gaussian noise generation
       void generateCircleData(
               const double cx,
               const double cy,
               const double r,
               const double noiseVar,
               std::vector<double> & x, std::vector<double> & y)
       {
           std::default_random_engine generator;
           std::normal_distribution<double> distribution(0,noiseVar);

           for(double i = 0 ; i < 3.1415*2.0 ; i += 3.1415/18.0){
               double xnoise = distribution(generator);
               double ynoise = distribution(generator);
               x.push_back(r*cos(double(i)) + cx); //cx
               y.push_back(r*sin(double(i)) + cy); //cy

               if (noiseVar != 0){
                   x[i] += xnoise;
                   y[i] += ynoise;
               }
           }
       }

       TEST(CircleFitting, idealCase)
       {
           std::vector<double> x;
           std::vector<double> y;

           double cx = 0; // circle center C:(cx,cy)
           double cy = 0;
           double radius = 1; //circle radius
           double noiseVar = 0.0;

           generateCircleData(cx,cy,radius,noiseVar,x,y);

           CircleFittingProblem * circleFitting = new CircleFittingProblem();
           circleFitting->setData(x,y);

           robest::LMedS solver;
           solver->solve(circleFitting);

           double res_cx,res_cy,res_r;
           circleFitting->getResult(res_cx,res_cy,res_r);

           ASSERT_NEAR(    cx, res_cx, 1.0e-11);
           ASSERT_NEAR(    cy, res_cy, 1.0e-11);
           ASSERT_NEAR(radius,  res_r, 1.0e-11);
       }   

       TEST(CircleFitting, idealCase2)
       {
           std::vector<double> x;
           std::vector<double> y;

           double cx = 24.8; // circle center C:(cx,cy)
           double cy = 8.10;
           double radius = 26.03; //circle radius
           double noiseVar = 0.0;

           generateCircleData(cx,cy,radius,noiseVar,x,y);

           CircleFittingProblem * circleFitting = new CircleFittingProblem();
           circleFitting->setData(x,y);

           double thres = 0.001;
           int nbIter = 20;
           robest::LMedS solver;
           solver->solve(circleFitting, thres, nbIter);

           double res_cx,res_cy,res_r;
           circleFitting->getResult(res_cx,res_cy,res_r);

           ASSERT_NEAR(    cx, res_cx, 1.0e-11);
           ASSERT_NEAR(    cy, res_cy, 1.0e-11);
           ASSERT_NEAR(radius,  res_r, 1.0e-11);
       }

       TEST(CircleFitting, smallNoise)
      {
           std::vector<double> x;
           std::vector<double> y;

           double cx = 3.552356; // circle center C:(cx,cy)
           double cy = 1.58452;
           double radius = 13.2548; //circle radius
           double noiseVar = 0.001;

           generateCircleData(cx,cy,radius,noiseVar,x,y);

           CircleFittingProblem * circleFitting = new CircleFittingProblem();
           circleFitting->setData(x,y);

           robest::LMedS solver;
           solver->solve(circleFitting);

           double res_cx,res_cy,res_r;
           circleFitting->getResult(res_cx,res_cy,res_r);

           ASSERT_NEAR(    cx, res_cx, 1.0e-3);
           ASSERT_NEAR(    cy, res_cy, 1.0e-3);
           ASSERT_NEAR(radius,  res_r, 1.0e-3);
       }

       TEST(CircleFitting, outliers)
       {
           std::vector<double> x = {1,0,-1, 0, sqrt(2)/2.0, 24,  8, 26};
           std::vector<double> y = {0,1, 0,-1, sqrt(2)/2.0,  8, 10,  3};

           CircleFittingProblem * circleFitting = new CircleFittingProblem();
           circleFitting->setData(x,y);

           robest::LMedS solver;
           auto nbIter = solver->calculateIterationsNb(x.size(),0.99,0.45);
           solver->solve(circleFitting, 0.1, nbIter);

           double res_cx,res_cy,res_r;
           circleFitting->getResult(res_cx,res_cy,res_r);

           ASSERT_NEAR( 0.0, res_cx, 1.0e-11);
           ASSERT_NEAR( 0.0, res_cy, 1.0e-11);
           ASSERT_NEAR( 1.0,  res_r, 1.0e-11);
       }

       TEST(CircleFitting, isDegenerate)
       {
           // Generate data
           // y = k*x + b
           std::vector<double> x1 = {0,1,2};
           std::vector<double> y1 = {0,1,2};

           // Define estimation problem
           CircleFittingProblem * circleFitting = new CircleFittingProblem();
           circleFitting->setData(x1,y1);
           ASSERT_TRUE(circleFitting->isDegenerate({0,1,2}));

           std::vector<double> x2 = {0,1,2};
           std::vector<double> y2 = {0,1.001,2};
           circleFitting->setData(x2,y2);
           ASSERT_TRUE(circleFitting->isDegenerate({0,1,2}));

           std::vector<double> x3 = {0,1,5};
           std::vector<double> y3 = {0,1,2};
           circleFitting->setData(x3,y3);
           ASSERT_TRUE(!circleFitting->isDegenerate({0,1,2}));
       }


Explanation
~~~~~~~~~~~

*Step 1: Declaration*
^^^^^^^^^^^^^^^^^^^^^

    In this step, we will demonstrate one of the ways to organize the structure of the header file for the circleFitting class.

    - Including libraries

    .. code-block:: c++

       #include <list>
       #include <vector>
       #include <cmath>
       #include <iostream>

       #include "robust_estim.hpp"

    - Defining global class parameters

    .. code-block:: c++

       // Creating a new data structure - 2D point
       struct Point2d{
           double x;
           double y;
       };
       
       // Definition of a new data type - vector of 2D points
       typedef std::vector<Point2d> Point2Dvector;

    - Inheritance from `EstimationProblem`

    .. code-block:: c++

       class CircleFittingProblem : public robest::EstimationProblem

    - Declaring public class attributes

    .. code-block:: c++

       public:
    	   // Class constructor
           CircleFittingProblem();

           // Class destructor
           ~CircleFittingProblem();

           // Data setting function
           void setData(std::vector<double> & x, std::vector<double> & y);

           // Residual calculation function
           double estimErrorForSample(int i);

           // Function of calculating the parameters of the model
           void   estimModelFromSamples(const std::vector<int> & samplesIdx);

           // Data size calculation function
           int getTotalNbSamples() const{
               return (int) points.size();
           }

           // Output function calculated parameters
           void getResult(double & res_cx, double & res_cy, double & res_r) const{
               res_cx = this->cx;
               res_cy = this->cy;
               res_r  = this->r;
           }

           // Function to check the correctness of the selected points
           bool isDegenerate(const std::vector<int> & samplesIdx);


    - Declaring private class attributes

    .. code-block:: c++

       private:

           Point2Dvector points; // input data
           double cx;	// x coordinate of the center	
           double cy;	// y coordinate of the center
           double r; 	// circle radius


*Step 2: Definition*
^^^^^^^^^^^^^^^^^^^^

    Once all the major class attributes are declared, the next step is to define them inside the source file.

    - Definition of the constructor and destructor

    .. code-block:: c++
       
       // Don't forget to include the header file!
       #include "CircleFitting.hpp"

       // Class constructor definition
       CircleFittingProblem::CircleFittingProblem(){

	   // Setting the number of parameters for the equation of a circle
           setNbParams(3);

	   // Setting the minimum number of points needed to calculate a circle model
           setNbMinSamples(3);
       }
       
       // Class destructor definition
       CircleFittingProblem::~CircleFittingProblem(){

       }

    - Function of setting the data

    .. code-block:: c++

       void CircleFittingProblem::setData(std::vector<double> & x, std::vector<double> & y)
       {
	   // Clearing the internal vector of 2D points
           points.clear();

           // Filling the internal vector of 2D points by external values 
           for (int i = 0; i < x.size(); i++){
               Point2d data;
               data.x=x[i];
               data.y=y[i];
               points.push_back(data);
           }
       }

    - Defining the function of calculating the residual

       The main task of this function is to calculate the distance between a given point and 
       the surface of a circle.This value is called the residual. Knowledge of this value is necessary
       for further error calculation using the loss function.

       What is the distance between a circle *C(x,y)* and and a point P(x\ :sub:`p`,y\ :sub:`p`) ?

       The equation of this circle is:

       .. math::
	  
	   (x - x_c)^2+(y - y_c)^2 = r^2

       where *x*\ :sub:`c` and *y*\ :sub:`c` are the coordinates of the centre of circle,
       and *r* is a radius of the circle.

       .. image:: images/distPointCircle.jpg
	   :width: 353px
	   :height: 358px
	   :scale: 75 %
	   :alt: calculated model
	   :align: center

       The ray *OP* , starting at the origin *O* and passing through the point *P*, 
       intersects the circle at the point closest to *P*. So, the distance between 
       the circle and the point will be the difference of the distance of the point 
       from the origin and the radius of the circle - *D*. 

       Using the Distance Formula between two point in Cartesian system of coordinates, 
       the final form of the residual function is:

       .. math::
	  
           D = |\sqrt{(x_p - x_c)^2 + (y_p - y_c)^2} - r|

       Thus, the function code for calculating the residual takes the following form:
       
       .. code-block:: c++

           inline double CircleFittingProblem::estimErrorForSample(int i)
           {
              // distance circle-point = abs(<distance point-center> - radius)
              Point2d & p = points[i];
              return std::abs(sqrt((p.x-cx)*(p.x-cx)+(p.y-cy)*(p.y-cy)) - r);
           }
       
    - Function to calculate the model

       The task of this function is to calculate the basic parameters of the circle model,
       which has the following form:
       
       .. math::
	  
	   (x_i - x_c)^2 + (y_i - y_c)^2 = r^2

       Here, the parameters are the coordinates of the center and the radius of a circle.

       It is easy to see that the equation has three unknown parameters, so in order to 
       determine these parameters, it is enough to know the coordinates of three points.
       
       After expanding and rearranging the terms, the new equation of a circle is expressed below. 
       
       .. math::
	  
	   x_i^2 + y_i^2 = 2x_ix_c + 2y_iy_c + r^2 - x_c^2 - y_c^2

       where *x*\ :sub:`i` and *y*\ :sub:`i` are the coordinates of *i*\ :sub:`th` point.

       Taking into account the fact that we need three points to determine the parameters of
       the model, this equation can now be expressed in vector/matrix notation:
           
       .. math::
           f = \begin{bmatrix}
                   x_1^2 + y_1^2 \\
                   x_2^2 + y_2^2 \\
                   x_3^2 + y_3^2 
               \end{bmatrix}
       
       .. math::
           A = \begin{bmatrix}
                   2x_1 & 2y_1 & 1 \\
                   2x_2 & 2y_2 & 1 \\
                   2x_3 & 2y_3 & 1 
               \end{bmatrix}

       .. math::
           p = \begin{bmatrix}
                   x_c \\
                   y_c \\
                   r^2 - x_c^2 - y_c^2
               \end{bmatrix}
       
       The *f*  vector, the *A* matrix, and the *p* vector represents the consolidated terms of 
       the expanded circle equation. 

       The new equation is seen below. We can calculate the circle’s parameters using the terms 
       in the *p*.

       .. math::
           
           f = Ap

           p = A^{-1}f

       
       Thus, the final form of the parameter calculation function will be as follows:

       .. code-block:: c++

          inline void CircleFittingProblem::estimModelFromSamples(const std::vector<int> & samplesIdx){
              
              // Validation of selected points
              if( !isDegenerate(samplesIdx)){

                  Point2d & P = points[samplesIdx[0]];
                  Point2d & V = points[samplesIdx[1]];
                  Point2d & K = points[samplesIdx[2]];

                  //calculation of the coefficients of the mediating straight lines
                  double a = -(V.x - P.x)/(V.y - P.y);
                  double b = (V.x * V.x - P.x * P.x + V.y * V.y - P.y * P.y)/(2* (V.y - P.y));
                  double c = -(K.x - V.x)/(K.y - V.y);
                  double d = (K.x * K.x - V.x * V.x + K.y * K.y - V.y * V.y)/(2* (K.y - V.y));

                  //calculate the coordinates of the center of the circle O(A,B)
                  cx = (b-d)/(c-a);
                  cy = a*cx + b;

                  //calculate the radius of a circle
                  r = sqrt((P.x - cx)*(P.x - cx)+(P.y - cy)*(P.y - cy));
              }
          }

    - Verification of a degenerate set of points

       Considering that the main estimation algorithm is based on the principle of 
       random selection of points, and for the correct construction of the model, 
       this set should consist of points that do not lie on one straight line. 
       We need to create a function that will check this condition.

       In order to check this condition, it is sufficient to determine whether 
       the two vectors, formed by these points, are collier or not.

       So, the function *isDegenerated* takes the following form:

       .. code-block:: c++

          inline bool CircleFittingProblem::isDegenerate(const std::vector<int> & samplesIdx)
          {
              Point2d & P = points[samplesIdx[0]];
              Point2d & V = points[samplesIdx[1]];
              Point2d & K = points[samplesIdx[2]];

              // verify that points P, V and K are not at the line -> verify that PV and PK are colinear:

  
              //1. calculate the directing coefficient of the line PV
              double f = (V.y-P.y)/(V.x-P.x);

              //2. calculate the directing coefficient of the line PK
              double h = (K.y-P.y)/(K.x-P.x);

              //3. PV and PK Are colineaire if and only if f = h
              return ( f - h < 1e-3 );
          }
	
*Step 3: Testing*
^^^^^^^^^^^^^^^^^

    Now that the `circleFitting` class is ready, you need to test it. In our library, we use Google tests.

    - Ideal cases
       
       The first two tests are aimed at general verification of the correctness of the class.
       We are testing the so-called ideal situation when the input set of points corresponds to
       a very simple example of a circle, there are no outliers and noises.

       .. code-block:: c++

          
          TEST(CircleFitting, idealCase)
          {
              std::vector<double> x;
              std::vector<double> y;

              double cx = 0; // circle center C:(cx,cy)
              double cy = 0;
              double radius = 1; //circle radius
              double noiseVar = 0.0;

              generateCircleData(cx,cy,radius,noiseVar,x,y);

              CircleFittingProblem * circleFitting = new CircleFittingProblem();
              circleFitting->setData(x,y);

              // Solving
              robest::LMedS solver;
              solver->solve(circleFitting);
	      
 	      // Getting the results
              double res_cx,res_cy,res_r;
              circleFitting->getResult(res_cx,res_cy,res_r);

              // Verifying the results
              ASSERT_NEAR(    cx, res_cx, 1.0e-11);
              ASSERT_NEAR(    cy, res_cy, 1.0e-11);
              ASSERT_NEAR(radius,  res_r, 1.0e-11);
          }   

          TEST(CircleFitting, idealCase2)
          {
              std::vector<double> x;
              std::vector<double> y;

              double cx = 24.8; // circle center C:(cx,cy)
              double cy = 8.10;
              double radius = 26.03; //circle radius
              double noiseVar = 0.0;

              generateCircleData(cx,cy,radius,noiseVar,x,y);

              CircleFittingProblem * circleFitting = new CircleFittingProblem();
              circleFitting->setData(x,y);

	      // Solving
              double thres = 0.001;
              int nbIter = 20;
              robest::LMedS solver;
              solver->solve(circleFitting, thres, nbIter);

	      // Getting the results
              double res_cx,res_cy,res_r;
              circleFitting->getResult(res_cx,res_cy,res_r);

	      // Verifying the results
              ASSERT_NEAR(    cx, res_cx, 1.0e-11);
              ASSERT_NEAR(    cy, res_cy, 1.0e-11);
              ASSERT_NEAR(radius,  res_r, 1.0e-11);
          }

    - Small noises case

       Now we will test the class on noisy data.

       .. code-block:: c++

          TEST(CircleFitting, smallNoise)
          {
              std::vector<double> x;
              std::vector<double> y;

              double cx = 3.552356; // circle center C:(cx,cy)
              double cy = 1.58452;
              double radius = 13.2548; //circle radius
              double noiseVar = 0.001;

              generateCircleData(cx,cy,radius,noiseVar,x,y);

              CircleFittingProblem * circleFitting = new CircleFittingProblem();
              circleFitting->setData(x,y);

	      // Solving
              robest::LMedS solver;
              solver->solve(circleFitting);

  	      // Getting the results
              double res_cx,res_cy,res_r;
              circleFitting->getResult(res_cx,res_cy,res_r);

	      // Verifying the results
              ASSERT_NEAR(    cx, res_cx, 1.0e-3);
              ASSERT_NEAR(    cy, res_cy, 1.0e-3);
              ASSERT_NEAR(radius,  res_r, 1.0e-3);
          }

    - Data with an outliers
       
       As the next text, we can check how effectively the algorithm copes with the corrupted data,
       in which there are outliers.

       .. code-block:: c++

          TEST(CircleFitting, outliers)
          {
              std::vector<double> x = {1,0,-1, 0, sqrt(2)/2.0, 24,  8, 26};
              std::vector<double> y = {0,1, 0,-1, sqrt(2)/2.0,  8, 10,  3};

              CircleFittingProblem * circleFitting = new CircleFittingProblem();
              circleFitting->setData(x,y);

	      // Solving
              robest::LMedS solver;
              auto nbIter = solver->calculateIterationsNb(x.size(),0.99,0.45);
              solver->solve(circleFitting, 0.1, nbIter);
	
	      // Getting the results
              double res_cx,res_cy,res_r;
              circleFitting->getResult(res_cx,res_cy,res_r);

              // Verifying the results
              ASSERT_NEAR( 0.0, res_cx, 1.0e-11);
              ASSERT_NEAR( 0.0, res_cy, 1.0e-11);
              ASSERT_NEAR( 1.0,  res_r, 1.0e-11);
          }

    - Degenerate case

       The last thing to check is the function of checking the correctness of the selected set of points.

       .. code-block:: c++

          TEST(CircleFitting, isDegenerate)
          {
              // Generate data
              // y = k*x + b
              std::vector<double> x1 = {0,1,2};
              std::vector<double> y1 = {0,1,2};

              // Define estimation problem
              CircleFittingProblem * circleFitting = new CircleFittingProblem();
              circleFitting->setData(x1,y1);

	      // Verifying a degenerate set 1
              ASSERT_TRUE(circleFitting->isDegenerate({0,1,2}));

              std::vector<double> x2 = {0,1,2};
              std::vector<double> y2 = {0,1.001,2};
              circleFitting->setData(x2,y2);

	      // Verifying a degenerate set 2
              ASSERT_TRUE(circleFitting->isDegenerate({0,1,2}));

              std::vector<double> x3 = {0,1,5};
              std::vector<double> y3 = {0,1,2};
              circleFitting->setData(x3,y3);

	      // Verifying a non-degenerate set
              ASSERT_TRUE(!circleFitting->isDegenerate({0,1,2}));
          }













