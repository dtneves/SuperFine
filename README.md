<h1>SuperFine</h1>
<!--------------------------------------------------------------------------------------------------------------------->
<h2>Fast and Accurate Supertree Estimation</h2>
<p>
    <b>SuperFine</b> is a meta-method that utilizes a novel two-step procedure in order to improve the accuracy and
    scalability of supertree methods. The paper <em>SuperFine: Fast and Accurate Supertree Estimation</em> (see [1])
    is, probably, the best introduction to <b>SuperFine</b>.
</p>
<p>
    SuperFine allows to create supertrees from source trees, as those provided in the <em>datasets</em> directory.
    There are two main methods, from which one can choose, to compute a supertree, those are:
    <ul>
        <li>Matrix Representation with Parsimony (<b>MRP</b>); and</li>
        <li>Matrix Representation with Likelihood (<b>MRL</b>).</li>
    </ul>
    To run <b>MRP</b> analyses, <em>PAUP*</em> has to be available
    on the system used to run the <code>runSuperFine.py</code> script. On the same way, to run <b>MRL</b> analyses,
    <em>FastTree</em> has also to be available on the system. In a Unix-like system, to find out if <em>PAUP*</em> and
    <em>FastTree</em> are available one may issue one of the following commands at a command line:
    <ul>
        <li><code>paup</code></li>
        <li><code>which paup</code></li>
        <li><code>FastTree</code></li>
        <li><code>which FastTree</code></li>
    </ul>
    Obviously, these binaries may exist but with different names (e.g. <em>PAUP</em>, <em>fasttree</em>, and so forth).
    Thus, is up to the user find these tools or, alternatively, make them available.
    However, one should notice that <b>SuperFine</b> expects <code>paup</code> and
    <code>FastTree</code> (case-sensitive!) binaries to be available.
</p>
<p>
    At the moment of this writing, 14.November.2016, this version of <b>SuperFine</b> was successfully tested with:
    <ul>
        <li>Python versions 2.7 and 3;</li>
        <li>PAUP* version 4a150; and</li>
        <li>FastTree version 2.1.9</li>
    </ul>
    To download <em>PAUP*</em> visit
    <a href="http://people.sc.fsu.edu/~dswofford/paup_test/">http://people.sc.fsu.edu/~dswofford/paup_test/</a>, and
    to download <em>FastTree</em> visit
    <a href="http://www.microbesonline.org/fasttree/">http://www.microbesonline.org/fasttree/</a>.
    Again, notice that <b>SuperFine</b> expects the binaries <em>paup</em> and
    <em>FastTree</em> to be available system-wide.
</p>
<p>
    If you use <b>SuperFine</b> in your research we would appreciate if you cite [1],
    you may also be interested in [2][3][4].
</p>
<!--------------------------------------------------------------------------------------------------------------------->
<h2>How to Run</h2>
<p>
    Passing the <code>-h</code> argument to the <code>runSuperFine.py</code> script allows to print
    a <em>help</em> message, which is useful to learn how to use <b>SuperFine</b>.
</p>
<p>
    <b>Usage Examples</b><br />
    <code>python ./runSuperFine.py -h</code><br />
    <code>python ./runSuperFine.py -r rmrp ./datasets/biological/seabirds/kennedy.source_trees_manual</code><br />
    <code>python ./runSuperFine.py -r rmrp ./datasets/simulated/100-taxa/50/sm_data.0.source_trees</code><br />
    <code>python ./runSuperFine.py -r fml ./datasets/biological/seabirds/kennedy.source_trees_manual</code><br />
    <code>python ./runSuperFine.py -r fml ./datasets/simulated/100-taxa/50/sm_data.0.source_trees</code>
</p>
<!--------------------------------------------------------------------------------------------------------------------->
<h2>References</h2>
<p>
    [1] M. S. Swenson, R. Suri, C. R. Linder, T. Warnow, SuperFine: Fast and Accurate Supertree Estimation,
        in: Systematic Biology, Vol. 61, Oxford University Press, 2012, pp. 214–227.<br />
    [2] D. T. Neves, T. Warnow, J. L. Sobral, K. Pingali, Parallelizing SuperFine,
        in: Proceedings of the 27th Annual ACM Symposium on Applied Computing, SAC ’12, ACM, 2012, pp. 1361–1367.<br />
    [3] D. T. Neves, J. L. Sobral, Towards a Faster and Accurate Supertree Inference,
        in: 20th IEEE Symposium on Computers and Communications, ISCC ’15, IEEE, 2015.<br />
    [4] D. T. Neves, J. L. Sobral, Parallel SuperFine - A Tool for Fast and Accurate Supertree Estimation:
        Features and Limitations, in: Future Generation Computer Systems, Elsevier, 2016.
</p>
