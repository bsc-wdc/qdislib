/* Add a small branded card to the bottom of Furo's right-hand column (below the
 * "On this page" table of contents): live GitHub / PyPI stats + quick links.
 * Injected in JS because Furo doesn't expose that slot through configuration. */
document.addEventListener("DOMContentLoaded", function () {
    var rail =
        document.querySelector(".toc-sticky") ||
        document.querySelector(".toc-drawer");
    if (!rail || rail.querySelector(".qd-rail-card")) {
        return;
    }
    var badge = function (alt, href, src) {
        return (
            '<a href="' +
            href +
            '" title="' +
            alt +
            '"><img loading="lazy" alt="' +
            alt +
            '" src="' +
            src +
            '"></a>'
        );
    };
    var card = document.createElement("div");
    card.className = "qd-rail-card";
    card.innerHTML = [
        '<div class="qd-rail-title">✂ Qdislib</div>',
        '<div class="qd-rail-badges">',
        badge(
            "GitHub stars",
            "https://github.com/bsc-wdc/qdislib",
            "https://img.shields.io/github/stars/bsc-wdc/qdislib?style=flat-square&logo=github&label=Stars&color=6f5aa6"
        ),
        badge(
            "GitHub forks",
            "https://github.com/bsc-wdc/qdislib/forks",
            "https://img.shields.io/github/forks/bsc-wdc/qdislib?style=flat-square&logo=github&label=Forks&color=13676b"
        ),
        badge(
            "PyPI downloads",
            "https://pepy.tech/project/qdislib",
            "https://img.shields.io/pepy/dt/qdislib?style=flat-square&logo=pypi&logoColor=white&label=Downloads&color=6f5aa6"
        ),
        "</div>",
        '<a href="https://doi.org/10.1145/3731599.3767547">📄 Read the paper (SC ’25)</a>',
        '<a href="https://bit.ly/bsc-wdc-community">💬 Join the community</a>',
    ].join("");
    rail.appendChild(card);
});
