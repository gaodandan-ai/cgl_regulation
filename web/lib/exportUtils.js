(function attachExportUtils(root, factory) {
    const moduleApi = factory(root);
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglExportUtils = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createExportUtils(root) {
    'use strict';

    function sanitizeSpreadsheetCell(value) {
        const text = value === null || value === undefined ? '' : String(value);
        return /^[=+\-@]/.test(text) ? `'${text}` : text;
    }

    function escapeCsvCell(value, { alwaysQuote = false, protectFormulas = true } = {}) {
        let text = value === null || value === undefined ? '' : String(value);
        if (protectFormulas) text = sanitizeSpreadsheetCell(text);
        const escaped = text.replace(/"/g, '""');
        return alwaysQuote || /[",\r\n]/.test(text) ? `"${escaped}"` : escaped;
    }

    function toCsv(rows, options = {}) {
        const lineEnding = options.lineEnding || '\n';
        const content = (rows || []).map(row =>
            (row || []).map(cell => escapeCsvCell(cell, options)).join(',')
        ).join(lineEnding);
        return options.bom ? `\uFEFF${content}` : content;
    }

    function download(content, filename, mime, options = {}) {
        const documentRef = options.documentRef || (root && root.document);
        const urlApi = options.urlApi || (root && root.URL);
        const BlobImpl = options.BlobImpl || (root && root.Blob);
        if (!documentRef || !urlApi || !BlobImpl) {
            throw new Error('Browser download APIs are unavailable');
        }
        const payload = options.bom ? ['\uFEFF', content] : [content];
        const blob = new BlobImpl(payload, { type: mime });
        const url = urlApi.createObjectURL(blob);
        const anchor = documentRef.createElement('a');
        anchor.href = url;
        anchor.download = filename;
        anchor.style.display = 'none';
        documentRef.body.appendChild(anchor);
        anchor.click();
        anchor.remove();
        setTimeout(() => urlApi.revokeObjectURL(url), 0);
    }

    return { sanitizeSpreadsheetCell, escapeCsvCell, toCsv, download };
});
