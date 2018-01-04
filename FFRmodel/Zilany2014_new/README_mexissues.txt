For macs, you only need xcode developer tools and not full Xcode.
However, you will run into a license acceptance issue. To fix that:
Edit the files under /Applications/MATLAB_RXXXX.app/bin/maci64/mexopts (probably there should be 3.xml file, all needs the same modification).

Locate the <XCODE_AGREED_VERSION> portion, comment the whole xml tag, e.g. wrap them with <!-- and --> like this:

<!--XCODE_AGREED_VERSION>
            <and diagnostic="Xcode is installed, but its license has not been accepted. Run Xcode and accept its license agreement." >
                <or>
                    <cmdReturns name="defaults read com.apple.dt.Xcode IDEXcodeVersionForAgreedToGMLicense"/>
                    <cmdReturns name="defaults read /Library/Preferences/com.apple.dt.Xcode IDEXcodeVersionForAgreedToGMLicense"/>
                </or>
                <cmdReturns name="&#10;agreed=$$ &#10; if echo $agreed | grep -E '[\.\&quot;]' >/dev/null; then &#10; lhs=`expr &quot;$agreed&quot; : '\([0-9]*\)[\.].*'` &#10;  rhs=`expr &quot;$agreed&quot; : '[0-9]*[\.]\(.*\)$'` &#10; if echo $rhs | grep -E '[\.&quot;]' >/dev/null; then &#10; rhs=`expr &quot;$rhs&quot; : '\([0-9]*\)[\.].*'` &#10; fi &#10; if [ $lhs -gt 4 ] || ( [ $lhs -eq 4 ] &amp;&amp; [ $rhs -ge 3 ] ); then &#10; echo $agreed &#10; else &#10; exit 1&#10; fi &#10; fi" />
            </and>
        </XCODE_AGREED_VERSION -->
Some notes:

These files are read only by default, you need to issue sudo chmod
644 * in that directory

Also, with the new clang, you will need to use 32-bit API for compiling, hence use
mex --compatibleArrayDims -v <C-files>
