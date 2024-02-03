<h2>1. Simulate ECD spectra</h2>
<p>Link：http://blog.molcalx.com.cn/2022/09/27/flare-ecd.html</p>

<h2>2. Torsion profiling</h2>
<p>Link:http://blog.molcalx.com.cn/2020/01/30/qm-torsion-profile.html</p>

<h2>3. Custom force field</h2>
<p>Link: http://blog.molcalx.com.cn/2020/12/17/custom-force-field.html</p>

<h2>4. Chemistry model</h2>
<p>Reference 1,2 and 3 </p>
<h2>5. Optimizaiton & Frequency calculation </h2>
<p>Optimization at BP86-D3/def2-tzvp level:</p>
<pre lang="python">
memory 24000 MB 

molecule triazole_zn {
2 1
    C            3.311141463382     1.388672716141    -0.132455063797
    N            1.884833345231     1.025177023894     0.043484026310
    C            1.254025091160    -0.054320546006    -0.457728777982
    N           -0.033045892061     0.033162706071    -0.033961339025
    N           -0.233960858151     1.139994388482     0.722214403874
    C            0.939263466825     1.740097045652     0.765668039406
    H            3.803704757611     0.634347732552    -0.772374626260
    H            3.368955725442     2.383933786760    -0.615462643757
    H            3.799389828311     1.413487992718     0.861690507697
    H            1.703430175308    -0.839624143660    -1.082877924916
    H            1.144432933625     2.680355668033     1.300493925735
    ZN          -1.619013907139    -0.992420975422    -0.232973752388
}

set {
    basis def2-tzvp
}

optimize('bp86-d3', dertype='energy')

</pre>
<h2>Reference</h2>
<ol>
  <li>PSI4 DFT functionals:<a href="https://psicode.org/psi4manual/4.0b5/dft_byfunctional.html">https://psicode.org/psi4manual/4.0b5/dft_byfunctional.html</a></li>
  <li>PSI4 Basis sets by elements: <a href="https://psicode.org/psi4manual/master/basissets_byelement.html">https://psicode.org/psi4manual/master/basissets_byelement.html</a></li>
  <li>M. Bursch, J.-M. Mewes, A. Hansen, S. Grimme, Angew. Chem. Int. Ed. 2022, 61, e202205735; Angew. Chem. 2022, 134, e202205735. <a href="https://onlinelibrary.wiley.com/doi/full/10.1002/anie.202205735">https://onlinelibrary.wiley.com/doi/full/10.1002/anie.202205735</a></li>
</ol>
