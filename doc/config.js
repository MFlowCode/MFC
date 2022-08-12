// This file is set as MATHJAX_CODEFILE in the Doxyfile. It configures how
// MathJax renders expressions in Markdown so that it is consistent with GitHub.

MathJax.Hub.Config({
    extensions: ["tex2jax.js"],
    jax: ["input/TeX", "output/HTML-CSS"],
    tex2jax: {
      inlineMath:  [ ['$',  '$'], ["\\(","\\)"] ],
      displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
      processEscapes: true
    },
    "HTML-CSS": {
      fonts: ["TeX"]
    }
});
