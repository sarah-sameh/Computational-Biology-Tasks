#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pyopenms import *
seq = AASequence.fromString("VAKA")
print("Total is =",seq.getMonoWeight())
sum=0
for aa in seq:
        sum +=aa.getMonoWeight()
        
print("Sum is =",sum)


# In[ ]:




