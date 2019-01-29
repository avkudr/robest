==================================
Introduction
==================================

At present, there is a large number of data estimation methods, but some of them
cannot guarantee a robustness of the obtained results. The main reason for inprecision is the presence of outliers or mismatches in the input data.
However, there is a family of algorithms working with such data: they are regrouped by the name of robust estimators
such as Least Median of Squared (LMedS) or RANSAC.

For example, below you can see the result of circle fitting using RANSAC applied to noisy data. Let's try to describe the algorithm step by step:

- *Step 1*: 
   
   The input of the algorithm is a quiet corrupted set of 2D points. The task is to define a circle model based on this set and avoid the effects of outliers.

   .. image:: images/exRANSACstep1.jpg
      :width: 556px
      :height: 494px
      :scale: 75 %
      :alt: step 1
      :align: center

- *Step 2*:

   In the next step, the algorithm estimates the model from a minimal needed data. 
   In case of circle, the minimum data needed to find its model (center coordinates and radius) is a set of three points. 
   Assume the algorithm chooses randomly three points (P\ :sub:`1`, P\ :sub:`2` and P\ :sub:`3`) and build a circle model on their basis.

   .. image:: images/exRANSACstep2.jpg
      :width: 556px
      :height: 494px
      :scale: 75 %
      :alt: step 2
      :align: center

- *Step 3*:

   The next step consists in calculation the loss function. In order words, we need to undestand how well this model fits all of the data.
   So, for each point in the set, an error is calculated: here we take the shortest distance between the point and the circle defined by the model. 
   If the error is below some predefined threshold value, then the point is classified as an inlier. Otherwise the point is considered an outlier.

   .. image:: images/exRANSACstep3.jpg
      :width: 556px
      :height: 494px
      :scale: 75 %
      :alt: step 3
      :align: center

- *The end*:

   Next, steps 2 and 3 are repeated multiple time in order to find the best model. For example, in case of RANSAC the best model 
   is the model that has the biggest number of inliers. 
   
   As a result, the algorithm will return the best selected model parameters. For a more detailed description of robust estimators 
   please refer to Algorithms section. 

   .. image:: images/exRANSACstep4.jpg
      :width: 556px
      :height: 494px
      :scale: 75 %
      :alt: step 4
      :align: center












