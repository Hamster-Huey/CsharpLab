# CsharpLab
## A look into good intentions with terrible practices 

I'm going to leave this code here just because it's my first ever major published project. 

I made this when I was a sophomore in college. My first major research job involved image processing and my professor wanted me to port matlab code to C# for use with a Kinect application (back when it v2 wasn't released yet, we had a dev kit). I wrote these functions to do image processing in C#, mainly writing just the exact functions I use in matlab for the project. They work, though I'm not sure how optimized my code is. It's cool to look back at however.

Look at how I just left the main file as Class1 
I didn't organize my functions based on... usage
I didn't need to make a class, could have used a namespace, as the functions don't rely on each other or any specific member variables. 
Lots of unnecessary stuff. 
