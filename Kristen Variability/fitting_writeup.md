
# Notes on fitting


## Accounting for size effects
Currently, it seems that the best fits that we get (from the perspective of AIC) are for size effects based on a Voigt profile, where the resonance position changes as a function of size, but the width (both $\gamma$ and $\sigma$) are not size dependent. This means we need to account for the size effects of the resonance position.

### starting

At present, the size effects are a shift from a starting resonance position, which is taken to be the bulk metal, in vacuum. We will call this $\mu_{metal}$.

### volume effects

Here the idea is that, as the particle shrinks, the g-value will change, and it will approach that of an isolated atom.  The deviation, to a first approximation might be expected to depend on the number of atoms present. When the number of atoms is infinite, then we expect to observe the resonance position associated with the metal.  When there is only 1 atom, we expect to observe the resonance position associated with a single atom.

We can express this as follows:

$\mu_{observed} = \mu_{metal} + \frac{\Delta\mu_{volume}}{N_{atoms}}$

The number of atoms present in the particle can be calculated by dividing the volume of the particle by the atomic volume of an atom:

$N_{atoms} = \frac{\frac{4}{3}\pi r_{particle}^3}{\frac{4}{3}\pi r_{atom}^3} = \frac{r_{particle}^3}{r_{atom}^3}$

### ligand effects
We also have ligands and solvent to worry about. Here, we might expect that we need to know how many ligands are present on the surface, and we might also want to consider that, as the number of atoms in the particle increases to infinity, the effect of the ligands will vanish. 

$\mu_{observed} = \mu_{metal} + \frac{\Delta\mu_{surface}\cdot f_{coverage}\cdot N_{surface}}{N_{atoms}}$

where $N_{surface}$ is the number of metal atoms at the surface and $f_{coverage}$ is the fraction of the metal atoms at the surface that have ligands bound to them.

We can estimate the number of atoms at the surface by calculating the volume of a shell that is the thickness of a single atom, and then dividing by the volume of a single atom. This is not perfect, as it is probably an over-estimate, but it is a place to start. 

$N_{surface} = \frac{\frac{4}{3}\pi [r_{particle}^3 - (r_{particle}-r_{atom})^3]}{\frac{4}{3}\pi r_{atom}^3}=\frac{r_{particle}^3 - (r_{particle}-2r_{atom})^3}{r_{atom}^3}$

so, then we have...

$\mu_{observed} = \mu_{metal} + \frac{\Delta\mu_{surface}\cdot f_{coverage}\cdot \frac{r_{particle}^3 - (r_{particle}-2r_{atom})^3}{r_{atom}^3}}{\frac{r_{particle}^3}{r_{atom}^3}} = \mu_{metal} + \frac{\Delta\mu_{surface}\cdot f_{coverage} [r_{particle}^3 - (r_{particle}-2r_{atom})^3]}{r_{particle}^3}$ 

## solvent effects
for the solvent effects, I again believe that it is largely dependent on the surface area to volume ratio. If we ascribe the solvent primarily to affecting the ligands, then we can just combine the solvent with the above, so that $\mu_{surface}$ is a combined ligand+solvent effect. If we think there is an additional effect, like what was described in our JACS paper, then we can try to include that as well. But for now, we have good fitting and so maybe we can stick with just combining them?


## End result

Combining the above thinking, we end up with the following equation for handing size effects:
$\mu_{observed} = \mu_{metal} + \mu_{volume}\frac{1}{N_{atoms}} + \mu_{surface}\frac{f_{coverage}\cdot N_{surface}}{N_{atoms}} = \mu_{metal} + \mu_{volume}\frac{r_{atom}^3}{r_{particle}^3} + \mu_{surface}\frac{f_{coverage}\cdot [r_{particle}^3 - (r_{particle} - 2r_{atom})^3]}{r_{particle}^3}$ 
Using this approach, we can fit a model, where each particle is described by a Lorentzian shape, where the resonance position of the Lorentzian changes according to the above, but the width of the Lorentzian does not.  This yields (in 30 seconds):
![[Pasted image 20250618122650.png]]

And if we then do the same, but shift to a Voigt profile, where the resonance position changes, but the $\gamma$ and $\sigma$ values do not, we obtain (in 40 minutes):

![[Pasted image 20250618152017.png]]
This is very good. I am not sure if the Gaussian component is 'real' or if it is just smoothing out the discrete counting we have, but at least there is a reasonable physical explanation for it.
# Next
- [ ] Compare this approach to the powder pattern fitting (lorentzian and voigt) 
- [ ] think about how to extract a meaningful value of resonance for a population?
	- do I use the mean value from a lognormal fit to talk about a mean value of the g-factor?
	- do I use the fit or do I try to calculate a mean value from the experimental data?  How to do this for lognormal?
- [ ] figure out what might be useful values for all these things.  Maybe convert the values we have into units that are physically relevant... like $g$ space and seconds. 
- [ ]  is it worth providing a guess for g_metal? Or even $\Delta g_{volume}$, if we know what the atom should be at?  Or is this really meaningful, as we get to a place that the electronic energy levels spread out too much? Maybe we no longer expect smooth transitions.