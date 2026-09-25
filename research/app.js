document.addEventListener('DOMContentLoaded', () => {

    // --- Research topic switcher (ARIA tabs with roving tabindex) ---
    const tabs = [...document.querySelectorAll('[role="tab"]')];
    const panes = document.querySelectorAll('.research-pane');
    const select = (tab, focus) => {
        tabs.forEach(t => {
            const on = t === tab;
            t.classList.toggle('active', on);
            t.setAttribute('aria-selected', String(on));
            t.tabIndex = on ? 0 : -1;
        });
        panes.forEach(p => p.classList.toggle('active', p.id === tab.dataset.target));
        if (focus) tab.focus();
        tab.scrollIntoView({ block: 'nearest', inline: 'center', behavior: 'smooth' });
    };
    tabs.forEach((tab, i) => {
        tab.addEventListener('click', () => select(tab, false));
        tab.addEventListener('keydown', e => {
            let next = null;
            if (e.key === 'ArrowRight') next = tabs[(i + 1) % tabs.length];
            else if (e.key === 'ArrowLeft') next = tabs[(i - 1 + tabs.length) % tabs.length];
            else if (e.key === 'Home') next = tabs[0];
            else if (e.key === 'End') next = tabs[tabs.length - 1];
            if (next) { e.preventDefault(); select(next, true); }
        });
    });

    // --- Drop the right-edge fade once the tab strip is scrolled to its end ---
    const wrap = document.querySelector('.segmented-wrap');
    if (wrap) {
        const onWrapScroll = () => wrap.classList.toggle('at-end', wrap.scrollLeft + wrap.clientWidth >= wrap.scrollWidth - 2);
        wrap.addEventListener('scroll', onWrapScroll, { passive: true });
        window.addEventListener('resize', onWrapScroll);
        onWrapScroll();
    }

    // --- Hairline under the nav once the page scrolls ---
    const navbar = document.querySelector('.navbar');
    const onScroll = () => navbar && navbar.classList.toggle('scrolled', window.scrollY > 8);
    window.addEventListener('scroll', onScroll, { passive: true });
    onScroll();

    // --- Highlight the nav link for the section in view ---
    const links = [...document.querySelectorAll('.nav-link[href^="#"]')];
    const sections = links.map(l => document.querySelector(l.getAttribute('href'))).filter(Boolean);
    if ('IntersectionObserver' in window && sections.length) {
        const io = new IntersectionObserver(entries => {
            entries.forEach(e => {
                if (!e.isIntersecting) return;
                links.forEach(l => l.classList.toggle('active', l.getAttribute('href') === '#' + e.target.id));
            });
        }, { rootMargin: '-45% 0px -50% 0px' });
        sections.forEach(s => io.observe(s));
    }
});
