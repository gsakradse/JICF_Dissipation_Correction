# RunDissipation_ArcLength_Package

Repository contains files to run the Structure function modification to the dissipation calculated with the gradient method.

# Potential issues:
    
* Normalization
    
    The script is set up to deal with the instantaneous velocity data with bad normalization, will have to adjust a normalization factor to deal with that if it is fixed  

* Vector spacing

    Within the ModifiedStructureFun file there is a variable called 'del', which is the vector spacing. It is probably different for the arc-length. We should discuss.