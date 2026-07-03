/* Add a small branded card to the bottom of Furo's right-hand column (below the
 * "On this page" TOC): live GitHub / PyPI stats + quick links. Injected in JS
 * because Furo doesn't expose that slot through configuration.
 *
 * Shields.io badges are static images, so they can't adapt to the page theme on
 * their own: we detect Furo's light/dark mode and request each badge with a
 * matching labelColor (white in light mode, near-black in dark mode), then
 * re-render when the theme toggle or the OS colour scheme changes. */
document.addEventListener("DOMContentLoaded", function () {
    var rail =
        document.querySelector(".toc-sticky") ||
        document.querySelector(".toc-drawer");
    if (!rail || rail.querySelector(".qd-rail-card")) {
        return;
    }

    var BADGES = [
        {
            alt: "GitHub stars",
            href: "https://github.com/bsc-wdc/qdislib",
            path: "github/stars/bsc-wdc/qdislib",
            label: "Stars",
            logo: "github",
            color: "6f5aa6",
        },
        {
            alt: "GitHub forks",
            href: "https://github.com/bsc-wdc/qdislib/forks",
            path: "github/forks/bsc-wdc/qdislib",
            label: "Forks",
            logo: "github",
            color: "13676b",
        },
        {
            alt: "PyPI downloads",
            href: "https://pepy.tech/project/qdislib",
            path: "pepy/dt/qdislib",
            label: "Downloads",
            logo: "pypi",
            color: "6f5aa6",
        },
    ];

    // Path to the docs root, derived from any _static asset URL, so the local
    // cheat-sheet link resolves on nested pages (notebooks/, api/...) too.
    var staticEl = document.querySelector(
        'link[href*="_static/"], script[src*="_static/"]'
    );
    var staticUrl = staticEl
        ? staticEl.getAttribute("href") || staticEl.getAttribute("src")
        : "";
    var root = staticUrl.slice(0, Math.max(0, staticUrl.indexOf("_static/")));

    var card = document.createElement("div");
    card.className = "qd-rail-card";
    card.innerHTML = [
        '<div class="qd-rail-title">✂ Qdislib</div>',
        '<div class="qd-rail-badges"></div>',
        '<a href="https://doi.org/10.1145/3731599.3767547">📄 Read the paper (SC ’25)</a>',
        '<a href="' + root + 'Qdislib_cheatsheet.pdf">📋 Cheat sheet (PDF)</a>',
        '<a href="https://bit.ly/bsc-wdc-community">💬 Join the community</a>',
    ].join("");
    rail.appendChild(card);
    var badgesEl = card.querySelector(".qd-rail-badges");

    function isDark() {
        var theme = document.body.dataset.theme; // "auto" | "light" | "dark"
        if (theme === "dark") return true;
        if (theme === "light") return false;
        return (
            window.matchMedia &&
            window.matchMedia("(prefers-color-scheme: dark)").matches
        );
    }

    function render() {
        var dark = isDark();
        var labelColor = dark ? "1b1b1b" : "ffffff"; // black in dark, white in light
        var logoColor = dark ? "ffffff" : "555555";
        badgesEl.innerHTML = BADGES.map(function (b) {
            var src =
                "https://img.shields.io/" +
                b.path +
                "?style=flat-square&logo=" +
                b.logo +
                "&logoColor=" +
                logoColor +
                "&label=" +
                b.label +
                "&labelColor=" +
                labelColor +
                "&color=" +
                b.color;
            return (
                '<a href="' +
                b.href +
                '" title="' +
                b.alt +
                '"><img loading="lazy" alt="' +
                b.alt +
                '" src="' +
                src +
                '"></a>'
            );
        }).join("");
    }

    render();

    // Re-render when Furo's theme toggle flips body[data-theme]...
    new MutationObserver(render).observe(document.body, {
        attributes: true,
        attributeFilter: ["data-theme"],
    });
    // ...and when the OS colour scheme changes while in "auto".
    if (window.matchMedia) {
        window
            .matchMedia("(prefers-color-scheme: dark)")
            .addEventListener("change", render);
    }
});
