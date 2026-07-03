/* Add a small branded "quick links" card to the bottom of Furo's right-hand
 * column (below the "On this page" table of contents). Injected in JS because
 * Furo doesn't expose that slot through configuration. */
document.addEventListener("DOMContentLoaded", function () {
    var rail = document.querySelector(".toc-sticky") || document.querySelector(".toc-drawer");
    if (!rail || rail.querySelector(".qd-rail-card")) {
        return;
    }
    var card = document.createElement("div");
    card.className = "qd-rail-card";
    card.innerHTML = [
        '<div class="qd-rail-title">✂ Qdislib</div>',
        '<a href="https://github.com/bsc-wdc/qdislib">⭐ Star on GitHub</a>',
        '<a href="https://doi.org/10.1145/3731599.3767547">📄 Read the paper (SC ’25)</a>',
        '<a href="https://qdislib.readthedocs.io">📘 Docs home</a>',
        '<a href="https://bit.ly/bsc-wdc-community">💬 Join the community</a>',
    ].join("");
    rail.appendChild(card);
});
