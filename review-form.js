// Leave-a-review form.
// Sends the review to Partho through FormSubmit (formsubmit.co), a free relay for static sites.
// First submission triggers a one-time activation email to the address below; click it once.
// If the relay fails (or is not yet activated), shows an "Email it instead" link and a copy button.
(() => {
  const form = document.getElementById('review-form');
  if (!form) return;
  const TO = 'pkpaul@uci.edu';
  const ENDPOINT = 'https://formsubmit.co/ajax/' + TO;
  const msg = document.getElementById('rf-msg');
  const copyBtn = document.getElementById('rf-copy');
  const mailLink = document.getElementById('rf-mail');
  const count = document.getElementById('rf-count');
  const btn = form.querySelector('button[type=submit]');
  const review = form.elements.review;
  review.addEventListener('input', () => { count.textContent = review.value.length; });
  form.addEventListener('input', () => { if (msg.classList.contains('err')) { msg.textContent = ''; msg.className = 'rf-msg'; } });

  const say = (text, kind) => { msg.textContent = text; msg.className = 'rf-msg' + (kind ? ' ' + kind : ''); };

  function fields() {
    const f = form.elements;
    return {
      course: f.subject.value,
      name: f.name.value.trim() || 'Anonymous',
      school: f.school.value.trim(),
      review: f.review.value.trim().replace(/\s*\n+\s*/g, ' '),
      consent: f.consent.checked ? 'Yes' : 'No'
    };
  }
  function offerFallback(d) {
    const subject = `Tutoring review · ${d.course}`;
    const body = [d.review, '', `Name to show: ${d.name}`, `Course or subject: ${d.course}`,
      d.school ? `School: ${d.school}` : null, `OK to post on the website: ${d.consent}`].filter(x => x !== null).join('\n');
    copyBtn.dataset.text = `To: ${TO}\nSubject: ${subject}\n\n${body}`;
    mailLink.href = `mailto:${TO}?subject=${encodeURIComponent(subject)}&body=${encodeURIComponent(body)}`;
    copyBtn.hidden = false; mailLink.hidden = false;
  }

  form.addEventListener('submit', async e => {
    e.preventDefault();
    if (form.elements._honey.value) return;               // bot
    if (!form.elements.subject.value) { say('Choose the course or subject first.', 'err'); form.elements.subject.focus(); return; }
    if (review.value.trim().length < 20) { say('Write at least a sentence or two.', 'err'); review.focus(); return; }
    const d = fields();
    btn.disabled = true; say('Sending…');
    try {
      const res = await fetch(ENDPOINT, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
        body: JSON.stringify({
          _subject: `Tutoring review · ${d.course}`,
          _template: 'table',
          _captcha: 'false',
          'Course or subject': d.course,
          'Name to show': d.name,
          'School': d.school || '—',
          'Review': d.review,
          'OK to post on website': d.consent
        })
      });
      const data = await res.json().catch(() => ({}));
      if (res.ok && String(data.success) === 'true') {
        form.reset(); count.textContent = '0';
        copyBtn.hidden = true; mailLink.hidden = true;
        say('Thank you. Your review was sent to Partho.', 'ok');
      } else {
        throw new Error(data.message || 'not delivered');
      }
    } catch (err) {
      say('Sorry, the review could not be sent right now. You can email it instead:', 'err');
      offerFallback(d);
    } finally {
      btn.disabled = false;
    }
  });

  copyBtn.addEventListener('click', async () => {
    try { await navigator.clipboard.writeText(copyBtn.dataset.text || ''); say(`Copied. Paste it into an email to ${TO}.`, 'ok'); }
    catch { say(`Copy failed. Please email ${TO} directly.`, 'err'); }
  });
})();
