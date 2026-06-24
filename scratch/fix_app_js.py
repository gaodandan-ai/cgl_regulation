import os

def main():
    filepath = os.path.join("web", "app.js")
    with open(filepath, "rb") as f:
        content_bytes = f.read()

    # Decode with ignore to strip any corrupted multi-byte sequences (errors='ignore')
    content_str = content_bytes.decode("utf-8", errors="ignore")

    # Define the new initAiSummaryFeature and triggerAiSummary block
    new_summary_features = """function initAiSummaryFeature() {
    const btnSaveKey = document.getElementById('btn-save-key');
    const btnClearKey = document.getElementById('btn-clear-key');
    const btnTriggerAi = document.getElementById('btn-trigger-ai');
    const apiKeyInput = document.getElementById('gemini-api-key-input');
    const keyConfigPanel = document.getElementById('ai-key-config-panel');
    const keyActivePanel = document.getElementById('ai-key-active-panel');
    
    // Multi-provider inputs
    const providerSelect = document.getElementById('ai-provider-select');
    const baseUrlInput = document.getElementById('ai-base-url-input');
    const modelInput = document.getElementById('ai-model-input');
    
    const customUrlWrapper = document.getElementById('ai-custom-url-wrapper');
    const modelWrapper = document.getElementById('ai-model-wrapper');
    const activeStatusText = document.getElementById('ai-active-status-text');

    const providerNames = {
        'google': 'Google Gemini',
        'openai': 'OpenAI',
        'deepseek': 'DeepSeek',
        'qwen': '通义千问',
        'kimi': 'Kimi',
        'zhipu': '智谱清言',
        'ollama': 'Ollama',
        'custom': '自定义接口'
    };

    const providerDefaults = {
        'google': { model: '', baseUrl: '' },
        'openai': { model: 'gpt-4o-mini', baseUrl: 'https://api.openai.com/v1' },
        'deepseek': { model: 'deepseek-chat', baseUrl: 'https://api.deepseek.com' },
        'qwen': { model: 'qwen-plus', baseUrl: 'https://dashscope.aliyuncs.com/compatible-mode/v1' },
        'kimi': { model: 'moonshot-v1-8k', baseUrl: 'https://api.moonshot.cn/v1' },
        'zhipu': { model: 'glm-4-flash', baseUrl: 'https://open.bigmodel.cn/api/paas/v4' },
        'ollama': { model: 'deepseek-r1', baseUrl: 'http://localhost:11434/v1' },
        'custom': { model: '', baseUrl: '' }
    };

    const hints = {
        'google': document.getElementById('ai-key-hint-google'),
        'openai': document.getElementById('ai-key-hint-openai'),
        'deepseek': document.getElementById('ai-key-hint-deepseek'),
        'qwen': document.getElementById('ai-key-hint-qwen'),
        'kimi': document.getElementById('ai-key-hint-kimi'),
        'zhipu': document.getElementById('ai-key-hint-zhipu'),
        'ollama': document.getElementById('ai-key-hint-ollama')
    };

    // Helper to toggle input fields visibility depending on selected provider
    function updateConfigFields() {
        const provider = providerSelect.value;
        
        // Hide all hints first
        Object.values(hints).forEach(h => {
            if (h) h.classList.add('hidden');
        });
        
        // Show current provider hint
        if (hints[provider]) {
            hints[provider].classList.remove('hidden');
        }

        // Toggle Base URL and Model visibility (hide only for Google Gemini)
        if (provider === 'google') {
            if (customUrlWrapper) customUrlWrapper.classList.add('hidden');
            if (modelWrapper) modelWrapper.classList.add('hidden');
        } else {
            if (customUrlWrapper) customUrlWrapper.classList.remove('hidden');
            if (modelWrapper) modelWrapper.classList.remove('hidden');
            
            // Adjust placeholders based on provider
            if (modelInput) {
                if (provider === 'custom') modelInput.placeholder = '例如: gpt-4o-mini';
                else modelInput.placeholder = `例如: ${providerDefaults[provider].model}`;
            }
        }

        // Adjust API Key label & requirements for Ollama
        const keyLabel = document.getElementById('ai-key-label');
        if (provider === 'ollama') {
            if (keyLabel) keyLabel.textContent = 'API 密钥 (Ollama 本地运行可选)';
            if (apiKeyInput) apiKeyInput.placeholder = '本地运行无需密钥，可为空...';
        } else {
            if (keyLabel) keyLabel.textContent = 'API 密钥 (API Key)';
            if (apiKeyInput) apiKeyInput.placeholder = '输入 API Key...';
        }
    }

    if (providerSelect) {
        providerSelect.addEventListener('change', () => {
            const provider = providerSelect.value;
            
            // Check if current inputs are empty or default values of ANY provider
            const currentModel = modelInput.value.trim();
            const currentBaseUrl = baseUrlInput.value.trim();
            
            const isModelDefaultOfAny = Object.values(providerDefaults).some(d => d.model === currentModel) || currentModel === '';
            const isBaseUrlDefaultOfAny = Object.values(providerDefaults).some(d => d.baseUrl === currentBaseUrl) || currentBaseUrl === '';
            
            if (isModelDefaultOfAny && providerDefaults[provider]) {
                modelInput.value = providerDefaults[provider].model;
            }
            if (isBaseUrlDefaultOfAny && providerDefaults[provider]) {
                baseUrlInput.value = providerDefaults[provider].baseUrl;
            }
            
            updateConfigFields();
        });
    }

    // 1. Migrate legacy key if present
    const legacyKey = localStorage.getItem('gemini_api_key');
    if (legacyKey && !localStorage.getItem('ai_api_key')) {
        localStorage.setItem('ai_api_key', legacyKey);
        localStorage.setItem('ai_provider', 'google');
        localStorage.removeItem('gemini_api_key'); // clear legacy
    }

    // 2. Load configurations on initialize
    function loadSavedConfig() {
        const savedKey = localStorage.getItem('ai_api_key');
        const savedProvider = localStorage.getItem('ai_provider') || 'google';
        const savedModel = localStorage.getItem('ai_model') || '';
        const savedBaseUrl = localStorage.getItem('ai_base_url') || '';

        if (providerSelect) providerSelect.value = savedProvider;
        if (modelInput) modelInput.value = savedModel;
        if (baseUrlInput) baseUrlInput.value = savedBaseUrl;
        
        updateConfigFields();

        if (savedKey || savedProvider === 'ollama') {
            keyConfigPanel.classList.add('hidden');
            keyActivePanel.classList.remove('hidden');
            if (activeStatusText) {
                const name = providerNames[savedProvider] || 'AI';
                activeStatusText.innerHTML = `<i class="fa-solid fa-circle-check"></i> ${name} 已就绪`;
            }
            btnTriggerAi.disabled = false;
        } else {
            keyConfigPanel.classList.remove('hidden');
            keyActivePanel.classList.add('hidden');
            btnTriggerAi.disabled = true;
        }
    }

    loadSavedConfig();

    // 3. Save settings listener
    btnSaveKey.addEventListener('click', () => {
        const provider = providerSelect.value;
        const key = apiKeyInput.value.trim();
        const model = modelInput.value.trim();
        const baseUrl = baseUrlInput.value.trim();

        if (!key && provider !== 'ollama') {
            alert('请输入 API 密钥！');
            return;
        }

        if (provider === 'custom' && !baseUrl) {
            alert('使用自定义服务商时，必须输入接口基址 (Base URL)！');
            return;
        }

        localStorage.setItem('ai_provider', provider);
        localStorage.setItem('ai_api_key', key);
        localStorage.setItem('ai_model', model);
        localStorage.setItem('ai_base_url', baseUrl);

        apiKeyInput.value = '';
        loadSavedConfig();
    });

    // 4. Clear config listener
    btnClearKey.addEventListener('click', () => {
        localStorage.removeItem('ai_api_key');
        localStorage.removeItem('ai_provider');
        localStorage.removeItem('ai_model');
        localStorage.removeItem('ai_base_url');

        // Reset input fields
        if (apiKeyInput) apiKeyInput.value = '';
        if (modelInput) modelInput.value = '';
        if (baseUrlInput) baseUrlInput.value = '';
        if (providerSelect) providerSelect.value = 'google';

        loadSavedConfig();
        
        const summaryCard = document.getElementById('ai-summary-result');
        if (summaryCard) {
            summaryCard.classList.add('hidden');
            summaryCard.innerHTML = '';
        }
    });

    // 5. AI Trigger listener
    btnTriggerAi.addEventListener('click', () => {
        triggerAiSummary();
    });
}

async function triggerAiSummary() {
    const btnTriggerAi = document.getElementById('btn-trigger-ai');
    const summaryCard = document.getElementById('ai-summary-result');
    
    const locus = document.getElementById('info-locus').textContent.trim();
    const name = document.getElementById('info-name').textContent.trim();
    const apiKey = localStorage.getItem('ai_api_key');
    const provider = localStorage.getItem('ai_provider') || 'google';
    const model = localStorage.getItem('ai_model') || '';
    const baseUrl = localStorage.getItem('ai_base_url') || '';
    
    if (!locus || locus === '-') {
        alert('请先选择一个基因进行分析。');
        return;
    }
    if (!apiKey && provider !== 'ollama') {
        alert('请先在上方配置您的 API Key。');
        return;
    }
    
    // Set loading state
    btnTriggerAi.disabled = true;
    summaryCard.classList.remove('hidden');
    summaryCard.classList.add('loading');
    summaryCard.innerHTML = `
        <div class="ai-spinner"></div>
        <span style="font-weight: 500;">正在检索 PubMed 文献并请求 AI 总结中...</span>
    `;
    
    try {
        const headers = {
            'X-AI-API-Key': apiKey || '',
            'X-AI-Provider': provider
        };
        if (model) headers['X-AI-Model'] = model;
        if (baseUrl) headers['X-AI-Base-URL'] = baseUrl;

        if (apiKey) {
            headers['X-Gemini-API-Key'] = apiKey;
        }

        const response = await fetch(`/api/summarize?gene=${locus}&name=${name}`, {
            headers: headers
        });
        
        if (!response.ok) {
            throw new Error(`HTTP 错误: ${response.status}`);
        }
        
        const result = await response.json();
        
        if (result.error) {
            throw new Error(result.error);
        }
        
        // Remove loading state
        summaryCard.classList.remove('loading');
        
        // Render summary text (with simple markdown parser)
        let htmlContent = parseMarkdownToHtml(result.summary);
        
        // Append papers if present
        if (result.papers && result.papers.length > 0) {
            htmlContent += `
                <div class="ai-sources-list">
                    <div class="ai-sources-title"><i class="fa-solid fa-book"></i> 参考 PubMed 文献 (\${result.papers.length} 篇)</div>
            `;
            
            result.papers.forEach(p => {
                htmlContent += `
                    <div class="ai-source-item">
                        <i class="fa-solid fa-file-lines"></i>
                        <a href="https://pubmed.ncbi.nlm.nih.gov/\${p.pmid}" target="_blank" class="ai-source-link" title="点击在 PubMed 查看原始文献">
                            \${p.title} (PMID: \${p.pmid}) <i class="fa-solid fa-arrow-up-right-from-square" style="font-size: 8px;"></i>
                        </a>
                    </div>
                `;
            });
            
            htmlContent += `</div>`;
        }
        
        summaryCard.innerHTML = htmlContent;
        
    } catch (err) {
        console.error(err);
        summaryCard.classList.remove('loading');
        summaryCard.innerHTML = `
            <div style="color: #ef4444; font-weight: 500; display: flex; align-items: center; gap: 6px; margin-bottom: 6px;">
                <i class="fa-solid fa-circle-exclamation"></i> 总结生成失败
            </div>
            <p style="font-size: 11px; color: var(--text-secondary); line-height: 1.4;">
                \${err.message || '未知网络错误，请检查您的 API Key 是否正确或网络连接状态。'}
            </p>
        `;
    } finally {
        btnTriggerAi.disabled = false;
    }
}"""

    # 1. First fix the syntax error: "input.value = difunction initAiSummaryFeature() {"
    # Should be replaced by closing the array iteration and starting the summary function
    syntax_old = "input.value = difunction initAiSummaryFeature() {"
    
    if syntax_old in content_str:
        print("Found syntax truncation error. Fixing it...")
        replacement = """input.value = displayLabel;
    });
}

function initAiSummaryFeature() {"""
        content_str = content_str.replace(syntax_old, replacement)

    # 2. Find the old initAiSummaryFeature + triggerAiSummary block and replace it.
    start_marker = "function initAiSummaryFeature() {"
    end_marker = "function parseMarkdownToHtml("
    
    start_idx = content_str.find(start_marker)
    end_idx = content_str.find(end_marker)
    
    if start_idx == -1 or end_idx == -1:
        print("ERROR: Could not find summary feature markers in app.js!")
        return
        
    content_str = content_str[:start_idx] + new_summary_features + "\n\n" + content_str[end_idx:]

    # 3. Replace validation in initAiPathwayFeature
    pathway_old = """        const apiKey = localStorage.getItem('ai_api_key') || localStorage.getItem('gemini_api_key');
        if (!apiKey) {
            alert('要使用 AI 通路分析，请先在右侧详情面板配置您的 API Key！');
            // Open the detail sidebar if it is collapsed, and highlight the key input
            toggleRightSidebar(true);
            const apiKeyInput = document.getElementById('gemini-api-key-input');
            if (apiKeyInput) {
                apiKeyInput.focus();
                apiKeyInput.style.border = '2px solid #6366f1';
                setTimeout(() => {
                    apiKeyInput.style.border = '1px solid var(--border-color)';
                }, 2000);
            }
            return;
        }

        const provider = localStorage.getItem('ai_provider') || 'google';
        const model = localStorage.getItem('ai_model') || '';
        const baseUrl = localStorage.getItem('ai_base_url') || '';"""
        
    pathway_new = """        const apiKey = localStorage.getItem('ai_api_key') || localStorage.getItem('gemini_api_key');
        const provider = localStorage.getItem('ai_provider') || 'google';
        const model = localStorage.getItem('ai_model') || '';
        const baseUrl = localStorage.getItem('ai_base_url') || '';

        if (!apiKey && provider !== 'ollama') {
            alert('要使用 AI 通路分析，请先在右侧详情面板配置您的 API Key！');
            // Open the detail sidebar if it is collapsed, and highlight the key input
            toggleRightSidebar(true);
            const apiKeyInput = document.getElementById('gemini-api-key-input');
            if (apiKeyInput) {
                apiKeyInput.focus();
                apiKeyInput.style.border = '2px solid #6366f1';
                setTimeout(() => {
                    apiKeyInput.style.border = '1px solid var(--border-color)';
                }, 2000);
            }
            return;
        }"""
        
    content_str = content_str.replace(pathway_old, pathway_new)

    # 4. Replace validation in initAiGeneFeature
    gene_old = """        const apiKey = localStorage.getItem('ai_api_key') || localStorage.getItem('gemini_api_key');
        if (!apiKey) {
            alert('要使用 AI 基因分析，请先在右侧详情面板配置您的 API Key！');
            // Open the detail sidebar if it is collapsed, and highlight the key input
            toggleRightSidebar(true);
            const apiKeyInput = document.getElementById('gemini-api-key-input');
            if (apiKeyInput) {
                apiKeyInput.focus();
                apiKeyInput.style.border = '2px solid #6366f1';
                setTimeout(() => {
                    apiKeyInput.style.border = '1px solid var(--border-color)';
                }, 2000);
            }
            return;
        }

        const provider = localStorage.getItem('ai_provider') || 'google';
        const model = localStorage.getItem('ai_model') || '';
        const baseUrl = localStorage.getItem('ai_base_url') || '';"""
        
    gene_new = """        const apiKey = localStorage.getItem('ai_api_key') || localStorage.getItem('gemini_api_key');
        const provider = localStorage.getItem('ai_provider') || 'google';
        const model = localStorage.getItem('ai_model') || '';
        const baseUrl = localStorage.getItem('ai_base_url') || '';

        if (!apiKey && provider !== 'ollama') {
            alert('要使用 AI 基因分析，请先在右侧详情面板配置您的 API Key！');
            // Open the detail sidebar if it is collapsed, and highlight the key input
            toggleRightSidebar(true);
            const apiKeyInput = document.getElementById('gemini-api-key-input');
            if (apiKeyInput) {
                apiKeyInput.focus();
                apiKeyInput.style.border = '2px solid #6366f1';
                setTimeout(() => {
                    apiKeyInput.style.border = '1px solid var(--border-color)';
                }, 2000);
            }
            return;
        }"""
        
    content_str = content_str.replace(gene_old, gene_new)

    # 5. Find and delete duplicate broken block at the end
    dup_pattern = "}网络"
    dup_idx = content_str.find(dup_pattern)
    resizer_idx = content_str.find("function initSidebarResizer()")
    
    if dup_idx != -1 and resizer_idx != -1:
        print(f"Removing duplicate block of length {resizer_idx - (dup_idx + 1)} bytes.")
        content_str = content_str[:dup_idx + 1] + "\n\n" + content_str[resizer_idx:]
    else:
        print("Warning: Could not find duplicate block patterns.")

    # Save as clean UTF-8
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content_str)
    print("app.js successfully updated and saved in UTF-8!")

if __name__ == "__main__":
    main()
