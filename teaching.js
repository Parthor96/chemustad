(() => {
  const resources = window.TEACHING_RESOURCES || [];
  const groups = document.querySelector('#resource-groups');
  const filters = document.querySelector('#course-filters');
  if (!groups || !filters) return;

  const SMS = 'sms:+19494782499?&body=';
  // Same test as .only-fine in style.css: mouse users get email, phones get a text.
  const FINE = window.matchMedia && matchMedia('(hover: hover) and (pointer: fine)').matches;
  const COURSE_NAMES = { I: 'Organic I', II: 'Organic II' };
  const KINDS = [
    ['problems', 'Worksheet'],
    ['solutions', 'Solutions'],
    ['videos', 'Video']
  ];

  function element(tag, text, className) {
    const node = document.createElement(tag);
    if (text) node.textContent = text;
    if (className) node.className = className;
    return node;
  }

  function slug(text) {
    return String(text).toLowerCase().replace(/&/g, ' ').replace(/[^a-z0-9]+/g, '-').replace(/^-+|-+$/g, '');
  }

  // Every topic keeps the number and anchor it has in resources.js, whatever the filter.
  const topics = resources.map((row, index) => {
    const links = [];
    for (const [key, fallback] of KINDS) {
      for (const item of row[key] || []) {
        if (!item || !item.url) continue;
        let url;
        try { url = new URL(item.url, document.baseURI); } catch (e) { continue; }
        if (!['https:', 'http:'].includes(url.protocol)) continue;
        const label = item.label || fallback;
        links.push({ label: label.toLowerCase().startsWith(fallback.toLowerCase()) ? label : `${fallback} · ${label}`, url });
      }
    }
    return { ...row, number: index + 1, id: `topic-${slug(row.topic)}`, links };
  });

  function linkNode(link, topic) {
    const a = element('a', link.label);
    a.href = link.url.href;
    a.setAttribute('aria-label', `${topic}: ${link.label}`);
    if (link.url.origin !== location.origin) { a.target = '_blank'; a.rel = 'noopener noreferrer'; }
    return a;
  }

  function topicRow(t) {
    const li = element('li', '', t.links.length ? 'topic-item has-links' : 'topic-item');
    li.id = t.id;
    const head = element('div', '', 'topic-head');
    head.append(element('span', String(t.number), 'topic-num'), element('span', t.topic, 'topic-name'));
    if (t.links.length) {
      li.append(head);
      const links = element('p', '', 'topic-links');
      t.links.forEach((link, i) => {
        if (i) links.append(element('span', ' · ', 'sep'));
        links.append(linkNode(link, t.topic));
      });
      li.append(links);
    } else {
      const ask = element('a', FINE ? 'Email ›' : 'Text ›', 'topic-ask');
      const topicText = t.topic.replace(/&/g, 'and');
      const body = `Hi Partho, I'd like help with ${topicText.charAt(0).toLowerCase() + topicText.slice(1)}.`;
      ask.href = FINE
        ? `mailto:pk@chemustad.org?subject=${encodeURIComponent('Tutoring: ' + t.topic)}&body=${encodeURIComponent(body)}`
        : SMS + encodeURIComponent(body).replace(/'/g, '%27');
      ask.setAttribute('aria-label', `${FINE ? 'Email' : 'Text'} Partho about ${t.topic}`);
      head.append(ask);
      li.append(head);
    }
    return li;
  }

  function render(course) {
    const selected = topics.filter(t => course === 'all' || t.course === course);
    groups.replaceChildren();
    for (const key of Object.keys(COURSE_NAMES)) {
      const list = selected.filter(t => t.course === key);
      if (!list.length) continue;
      const section = element('section', '', 'topic-group');
      const h = element('h3', COURSE_NAMES[key]);
      h.id = `course-${key}-heading`;
      section.setAttribute('aria-labelledby', h.id);
      const ol = element('ol', '', 'topic-list');
      list.forEach(t => ol.append(topicRow(t)));
      section.append(h, ol);
      groups.append(section);
    }
    for (const button of filters.querySelectorAll('button')) {
      button.setAttribute('aria-pressed', String(button.dataset.course === course));
    }
  }

  // Honest status line, computed from the data.
  const sub = document.querySelector('#resources-sub');
  const posted = topics.filter(t => t.links.length).length;
  if (sub && topics.length) {
    sub.textContent = `${posted} of ${topics.length} topics ${posted === 1 ? 'has' : 'have'} material so far; more are posted as they are finished. ` +
      `Use ${FINE ? 'Email' : 'Text'} on any topic to ask Partho about it. Practice material is organic chemistry only for now.`;
  }

  // The course filter only earns its place once several topics have material.
  filters.hidden = topics.filter(t => t.links.length).length < 4;
  groups.hidden = false;
  filters.addEventListener('click', event => {
    const button = event.target.closest('button[data-course]');
    if (button) render(button.dataset.course);
  });
  render('all');

  // Topics are rendered by script, so jump to a #topic-… link once they exist.
  if (location.hash.startsWith('#topic-')) {
    const target = document.getElementById(decodeURIComponent(location.hash.slice(1)));
    if (target) target.scrollIntoView();
  }
})();
