/**
 * onboardingTour.js
 * Interactive Onboarding Guided Tour Engine for Cgl Regulation Explorer
 */
(function () {
    'use strict';

    const tourSteps = [
        {
            target: '.gene-input-wrapper',
            title: '1. 基因搜索与快捷示例',
            content: '在搜索框中输入 Corynebacterium glutamicum 的基因 locus tag (如 cg0350, cg0444) 或通用名称 (如 sigH, whiB4)，点击 Analyze 即可开启调控网络！',
            position: 'right'
        },
        {
            target: '.checkbox-group',
            title: '2. 调控边与过滤控制',
            content: '您可以勾选或取消勾选不同的边类型（激活 +、抑制 -、双重调控、sRNA 调控、PPI 蛋白质相互作用），实时过滤网络视图。',
            position: 'right'
        },
        {
            target: '#workflow-entrybar',
            title: '3. 顶栏 4 大工作流导航',
            content: '快速在【基因网络 Network】、【iModulon 模块】、【工程改造靶点 Priority】与【代谢途径 Pathways】四大核心工作流之间无缝切换。',
            position: 'bottom'
        },
        {
            target: '#right-sidebar',
            title: '4. AI 智能基因助手与文献挖掘',
            content: '右侧侧边栏集成了基于 RAG 大模型与 PubMed/Europe PMC 的智能基因问答助手，深度分析基因功能与参考文献。',
            position: 'left'
        },
        {
            target: '.quick-examples',
            title: '5. 快速上车，开启首次探索！',
            content: '点击下方任意示例基因标签（如 sigH / cg0350），立即可生成全景基因调控网络图谱！',
            position: 'bottom'
        }
    ];

    let currentStepIdx = 0;
    let isActive = false;
    let backdropEl = null;
    let popoverEl = null;

    function initElements() {
        if (backdropEl && popoverEl) return;

        backdropEl = document.createElement('div');
        backdropEl.id = 'tour-backdrop';
        backdropEl.className = 'tour-spotlight-backdrop';

        popoverEl = document.createElement('div');
        popoverEl.id = 'tour-popover';
        popoverEl.className = 'tour-popover';

        document.body.appendChild(backdropEl);
        document.body.appendChild(popoverEl);

        backdropEl.addEventListener('click', (e) => {
            if (e.target === backdropEl) {
                stopTour();
            }
        });
    }

    function positionPopover(targetEl, position) {
        if (!popoverEl) return;

        if (!targetEl) {
            // Fallback to screen center if target element is not found
            const viewportWidth = window.innerWidth;
            const viewportHeight = window.innerHeight;
            popoverEl.style.top = `${Math.max(20, (viewportHeight / 2) - 100)}px`;
            popoverEl.style.left = `${Math.max(20, (viewportWidth / 2) - 170)}px`;
            return;
        }

        const rect = targetEl.getBoundingClientRect();
        const popoverRect = popoverEl.getBoundingClientRect();
        const viewportWidth = window.innerWidth;
        const viewportHeight = window.innerHeight;

        let top = 0;
        let left = 0;

        switch (position) {
            case 'right':
                top = rect.top + (rect.height / 2) - (popoverRect.height / 2);
                left = rect.right + 18;
                if (left + popoverRect.width > viewportWidth - 20) {
                    left = rect.left - popoverRect.width - 18;
                }
                break;
            case 'left':
                top = rect.top + (rect.height / 2) - (popoverRect.height / 2);
                left = rect.left - popoverRect.width - 18;
                if (left < 20) {
                    left = rect.right + 18;
                }
                break;
            case 'top':
                top = rect.top - popoverRect.height - 18;
                left = rect.left + (rect.width / 2) - (popoverRect.width / 2);
                break;
            case 'bottom':
            default:
                top = rect.bottom + 18;
                left = rect.left + (rect.width / 2) - (popoverRect.width / 2);
                if (top + popoverRect.height > viewportHeight - 20) {
                    top = rect.top - popoverRect.height - 18;
                }
                break;
        }

        // Clamp inside screen viewport
        top = Math.max(20, Math.min(viewportHeight - popoverRect.height - 20, top));
        left = Math.max(20, Math.min(viewportWidth - popoverRect.width - 20, left));

        popoverEl.style.top = `${top}px`;
        popoverEl.style.left = `${left}px`;
    }

    function renderStep(idx) {
        if (idx < 0 || idx >= tourSteps.length) {
            stopTour();
            return;
        }

        // Auto-expand left sidebar if collapsed when targeting sidebar elements
        if (idx === 0 || idx === 1) {
            const sidebar = document.getElementById('sidebar');
            if (sidebar && sidebar.classList.contains('collapsed')) {
                sidebar.classList.remove('collapsed');
            }
        }

        currentStepIdx = idx;
        const step = tourSteps[idx];
        const targetEl = document.querySelector(step.target);

        // Remove active pulsing highlight from previous elements
        document.querySelectorAll('.tour-highlight-active').forEach(el => el.classList.remove('tour-highlight-active'));

        // Highlight spotlight position
        if (targetEl) {
            targetEl.classList.add('tour-highlight-active');
            try {
                targetEl.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
            } catch (e) {}

            const rect = targetEl.getBoundingClientRect();
            backdropEl.style.clipPath = `polygon(
                0% 0%, 0% 100%,
                ${rect.left - 10}px 100%,
                ${rect.left - 10}px ${rect.top - 10}px,
                ${rect.right + 10}px ${rect.top - 10}px,
                ${rect.right + 10}px ${rect.bottom + 10}px,
                ${rect.left - 10}px ${rect.bottom + 10}px,
                ${rect.left - 10}px 100%,
                100% 100%, 100% 0%
            )`;
        } else {
            backdropEl.style.clipPath = 'none';
        }

        const isFirst = idx === 0;
        const isLast = idx === tourSteps.length - 1;

        popoverEl.innerHTML = `
            <div class="tour-header">
                <span class="tour-step-badge">Step ${idx + 1} of ${tourSteps.length}</span>
                <button class="tour-close-btn" id="tour-close-action" title="Skip tour"><i class="fa-solid fa-xmark"></i></button>
            </div>
            <h4 class="tour-title">${step.title}</h4>
            <p class="tour-content">${step.content}</p>
            <div class="tour-footer">
                <button class="tour-btn tour-btn-secondary" id="tour-prev-action" ${isFirst ? 'disabled' : ''}>
                    <i class="fa-solid fa-chevron-left"></i> 上一步
                </button>
                <div class="tour-footer-right">
                    <button class="tour-btn tour-btn-link" id="tour-skip-action">跳过导览</button>
                    <button class="tour-btn tour-btn-primary" id="tour-next-action">
                        ${isLast ? '开始体验 <i class="fa-solid fa-rocket"></i>' : '下一步 <i class="fa-solid fa-chevron-right"></i>'}
                    </button>
                </div>
            </div>
        `;

        popoverEl.classList.add('show');
        backdropEl.classList.add('show');

        // Re-position after render
        requestAnimationFrame(() => {
            positionPopover(targetEl, step.position);
        });

        // Bind internal buttons
        document.getElementById('tour-prev-action')?.addEventListener('click', () => renderStep(currentStepIdx - 1));
        document.getElementById('tour-next-action')?.addEventListener('click', () => {
            if (isLast) {
                stopTour();
                // Trigger quick example search for instant AHA moment!
                if (typeof window.querySingleGene === 'function') {
                    window.querySingleGene('sigH');
                }
            } else {
                renderStep(currentStepIdx + 1);
            }
        });
        document.getElementById('tour-skip-action')?.addEventListener('click', stopTour);
        document.getElementById('tour-close-action')?.addEventListener('click', stopTour);
    }

    function startTour() {
        initElements();
        isActive = true;
        renderStep(0);
    }

    function stopTour() {
        isActive = false;
        document.querySelectorAll('.tour-highlight-active').forEach(el => el.classList.remove('tour-highlight-active'));
        if (backdropEl) backdropEl.classList.remove('show');
        if (popoverEl) popoverEl.classList.remove('show');
    }

    // Auto-bind click listeners when DOM is ready
    document.addEventListener('DOMContentLoaded', () => {
        const tourBtn = document.getElementById('btn-start-tour');
        if (tourBtn) {
            tourBtn.addEventListener('click', (e) => {
                e.preventDefault();
                e.stopPropagation();
                startTour();
            });
        }
    });

    // Expose API
    window.OnboardingTour = {
        start: startTour,
        stop: stopTour,
        next: () => renderStep(currentStepIdx + 1),
        prev: () => renderStep(currentStepIdx - 1)
    };
})();
