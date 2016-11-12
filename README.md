<h1>SuperFine</h1>
<!--------------------------------------------------------------------------------------------------------------------->
<h2>Fast and Accurate Supertree Estimation</h2>
<p>
    <b>SuperFine</b> is a meta-method that utilizes a novel two-step procedure in order to improve the accuracy and
    scalability of supertree methods. The paper <em>SuperFine: Fast and Accurate Supertree Estimation</em> (see [1])
    is, probably, the best introduction to <b>SuperFine</b>.
</p>
<p>
    To create a supertree from a set of source trees (as those provided in the <em>datasets</em> directory)
    using the Matrix Representation with Parsimony (MRP) supertree method requires <em>PAUP*</em> to be available
    on the system used to run the <code>runSuperFine.py</code> script. In a Unix-like system, to find out
    if <em>PAUP*</em> is available one may issue one of the following commands at a command line:<br />
    <code>paup</code><br />
    <code>which paup</code>
</p>
<p>
    This version of <b>SuperFine</b> was successfully tested with version <em>4a150</em> of <em>PAUP*</em>,
    which is, at the moment of this writing (18.October.2016), freely available at
    <a href=http://people.sc.fsu.edu/~dswofford/paup_test/>http://people.sc.fsu.edu/~dswofford/paup_test/</a>.
    Notice that <b>SuperFine</b> expects a binary named <em>paup</em> to be available system-wide.
</p>
<p>
    If you use SuperFine in your research we would appreciate if you cite [1],
    you may also be interested in [2], [3], and [4].
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
    <code>python ./runSuperFine.py -r rmrp ./datasets/simulated/100-taxa/50/sm_data.0.source_trees</code>
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
        Features and Limitations,
        in: Future Generation Computer Systems, Elsevier, 2016.
</p>
