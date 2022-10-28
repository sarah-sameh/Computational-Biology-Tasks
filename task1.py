#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('pip install pyopenms')


# In[5]:


import pyopenms


# In[6]:


print ("Avogadro's number is", pyopenms.Constants.AVOGADRO)


# In[7]:


from pyopenms import *

edb = ElementDB()

edb.hasElement("O")


# In[5]:


carbon = edb.getElement("C")


print(carbon.getName())

print(carbon.getSymbol())

print(carbon.getMonoWeight())

print(carbon.getAverageWeight())

print ("One mole of carbon weighs", 2*carbon.getAverageWeight(), "grams")


# In[8]:


hydrogen = edb.getElement("H")


print(hydrogen.getName())



print(hydrogen.getSymbol())



print(hydrogen.getMonoWeight())



print(hydrogen.getAverageWeight())



isotopes = hydrogen.getIsotopeDistribution()



print ("One mole of 16O2 weighs", 2*hydrogen.getMonoWeight(), "grams")


# In[9]:


edb = ElementDB()


oxygen_isoDist = {"mass": [], "abundance": []}

oxygen = edb.getElement("O")


isotopes = oxygen.getIsotopeDistribution()


for iso in isotopes.getContainer():
    print ("Oxygen isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    oxygen_isoDist["mass"].append(iso.getMZ())
    oxygen_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[10]:


sulfur_isoDist = {"mass": [], "abundance": []}


sulfur = edb.getElement("S")


isotopes = sulfur.getIsotopeDistribution()


for iso in isotopes.getContainer():
    print ("Sulfur isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    sulfur_isoDist["mass"].append(iso.getMZ())
    sulfur_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[11]:


import math
from matplotlib import pyplot as plt


def adjustText(x1, y1, x2, y2):
    if y1 > y2:
        plt.annotate('%0.3f' % (y2), xy=(x2, y2), xytext=(x2+0.5,y2+9),
                     textcoords='data',
                     arrowprops=dict(arrowstyle="->", color='r', lw=0.5),
                     horizontalalignment='right', verticalalignment='top')
    else:
        plt.annotate('%0.3f' % (y1), xy=(x1, y1), xytext=(x1+0.5,y1+9),
                     textcoords='data',
                     arrowprops=dict(arrowstyle="->", color='r', lw=0.5),
                     horizontalalignment='right', verticalalignment='top')


def plotDistribution(distribution):
    n = len(distribution["mass"])
    for i in range(0, n):
        plt.vlines(x=distribution["mass"][i], ymin=0, ymax=distribution["abundance"][i])
        if int(distribution["mass"][i - 1]) == int(distribution["mass"][i])                 and i != 0:
            adjustText(distribution["mass"][i - 1], distribution["abundance"][i - 1],
                       distribution["mass"][i], distribution["abundance"][i])
        else:
            plt.text(x=distribution["mass"][i],
                     y=(distribution["abundance"][i] + 2),
                     s='%0.3f' % (distribution["abundance"][i]), va='center',
                     ha='center')
    plt.ylim([0, 110])
    plt.xticks(range(math.ceil(distribution["mass"][0]) - 2,
                     math.ceil(distribution["mass"][-1]) + 2))


plt.figure(figsize=(10,7))

plt.subplot(1,2,1)
plt.title("Isotopic distribution of oxygen")
plotDistribution(oxygen_isoDist)
plt.xlabel("Atomic mass (u)")
plt.ylabel("Relative abundance (%)")

plt.subplot(1,2,2)
plt.title("Isotopic distribution of sulfur")
plotDistribution(sulfur_isoDist)
plt.xlabel("Atomic mass (u)")
plt.ylabel("Relative abundance (%)")

plt.show()


# In[12]:


edb = ElementDB()


isotopes = edb.getElement("C").getIsotopeDistribution().getContainer()


carbon_isotope_difference = isotopes[1].getMZ() - isotopes[0].getMZ()


isotopes = edb.getElement("N").getIsotopeDistribution().getContainer()



nitrogen_isotope_difference = isotopes[1].getMZ() - isotopes[0].getMZ()


# In[13]:


print ("Mass difference between 12C and 13C:", carbon_isotope_difference)


# In[14]:


print ("Mass difference between 14N and N15:", nitrogen_isotope_difference)


# In[15]:


print ("Relative deviation:", 100*(carbon_isotope_difference -nitrogen_isotope_difference)/carbon_isotope_difference, "%")


# In[16]:


methanol = EmpiricalFormula("CH3OH")


water = EmpiricalFormula("H2O")



ethanol = EmpiricalFormula("CH2") + methanol



print("Ethanol chemical formula:", ethanol.toString())


# In[17]:


print("Ethanol composition:", ethanol.getElementalComposition())


# In[18]:


print("Ethanol has", ethanol.getElementalComposition()[b"H"], "hydrogen atoms")


# In[19]:


lys = ResidueDB().getResidue("Lysine")


# In[20]:



print(lys.getName())


print(lys.getThreeLetterCode())


print(lys.getOneLetterCode())


print(lys.getAverageWeight())


print(lys.getMonoWeight())


print(lys.getPka())


print(lys.getFormula().toString())


# In[21]:


ox = ModificationsDB().getModification("Oxidation")


print(ox.getUniModAccession())


print(ox.getUniModRecordId())


print(ox.getDiffMonoMass())


print(ox.getId())


print(ox.getFullId())


print(ox.getFullName())


print(ox.getDiffFormula())


# In[24]:


isotopes = ox.getDiffFormula().getIsotopeDistribution(CoarseIsotopePatternGenerator(5))
for iso in isotopes.getContainer():
    print (iso.getMZ(), ":", iso.getIntensity())


# In[25]:


uridine = RibonucleotideDB().getRibonucleotide(b"U")



print(uridine.getName())


print(uridine.getCode())

print(uridine.getAvgMass())

print(uridine.getMonoMass())


print(uridine.getFormula().toString())

print(uridine.isModified())

methyladenosine = RibonucleotideDB().getRibonucleotide(b"m1A")


print(methyladenosine.getName())

print(methyladenosine.isModified())


# In[ ]:




