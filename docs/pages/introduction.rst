==================================
Introduction
==================================

At present, there are a large number of data estimation methods, but some of them
cannot guarantee a robustness of the results obtained. The main reason is the presence of outliers or mismatches in the data.
Despite this, among all these methods we can distinguish a family of solutions regrouped by the name of robust estimators
such as Least Median of Squared (LMedS) or RANSAC.

Consider the work of the method RANSAC on the example of the circle fitting.

- *Step I*: 
   
   The input of the algorithm is a quiet corrupted set of 2D points. The task is to define a circle model based on this set and avoid the effects of outliers.

   .. image:: images/exRANSACstep1.jpg
      :width: 556px
      :height: 494px
      :scale: 75 %
      :alt: step 1
      :align: center

- *Step II*:

   In the next step, the algorithm estimate the data. He choose randomly three points (P\ :sub:`1`, P\ :sub:`2` and P\ :sub:`3`) and build a circle model on their basis.

   .. image:: images/exRANSACstep2.jpg
      :width: 556px
      :height: 494px
      :scale: 75 %
      :alt: step 2
      :align: center

- *Step III*:

   The next step in the algorithm is to calculate the loss function based on a predetermined threshold value.
   First, for each point in the set, an error is calculated; this error is the calculation of the shortest distance
   between the point and the surface of the constructed circle model. If the square of the error is less than or 
   equal to the square of the threshold value, then the point is classified as an inlier, if not then the point is considered an outlier.

   The best estimation is the one at which the number of inliers is maximum.

   .. image:: images/exRANSACstep3.jpg
      :width: 556px
      :height: 494px
      :scale: 75 %
      :alt: step 3
      :align: center

- *Step IV*:

   All steps of the algorithm are repeated until the number of iterations reaches the specified value. 
   As a result, the algorithm will return the best selected model parameters.

   .. image:: images/exRANSACstep4.jpg
      :width: 556px
      :height: 494px
      :scale: 75 %
      :alt: step 4
      :align: center


   
      :align: center








