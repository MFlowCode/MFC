// This file is set as MATHJAX_CODEFILE in the Doxyfile. It configures how
// MathJax renders expressions in Markdown so that it is consistent with GitHub.

MathJax.Hub.Config({
    extensions: ["tex2jax.js"],
    jax: ["input/TeX", "output/HTML-CSS"],
    tex2jax: {
      inlineMath:  [ ['$',  '$'], ["\\(","\\)"] ],
      displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
      processEscapes: true,
      ignoreClass: "line" // Ignore code blocks: https://web.archive.org/web/20120430100225/http://www.mathjax.org/docs/1.1/options/tex2jax.html
    },
    "HTML-CSS": {
      fonts: ["TeX"]
    }
});
