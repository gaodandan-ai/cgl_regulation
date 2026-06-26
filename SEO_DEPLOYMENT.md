# SEO Deployment Checklist

This project is SEO-ready at the static page level, but search engines can only index it after it is deployed to a public HTTPS URL.

## 1. Set The Public URL

The current SEO files use:

```text
https://cgl-regulation.vercel.app/
```

If you deploy to another domain, replace this URL in:

- `web/index.html`
- `web/robots.txt`
- `web/sitemap.xml`

## 2. Deploy Publicly

Recommended options:

- Vercel, using the existing `vercel.json`
- Render/Fly/Railway/Docker, using `Dockerfile`
- GitHub Pages only if API routes are not needed

The site must be reachable by search engines over HTTPS. Localhost URLs such as `http://127.0.0.1:8000/index.html` cannot be indexed.

## 3. Submit For Indexing

After deployment:

1. Open Google Search Console.
2. Add the deployed domain or URL prefix.
3. Submit `https://your-domain/sitemap.xml`.
4. Use URL Inspection on the homepage and request indexing.
5. Repeat in Bing Webmaster Tools.

## 4. Verify

Check these public URLs:

```text
https://your-domain/
https://your-domain/robots.txt
https://your-domain/sitemap.xml
```

The homepage should include title, description, canonical URL, Open Graph tags, and JSON-LD structured data.
