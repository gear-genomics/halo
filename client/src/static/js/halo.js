import axios from 'axios'
import pako from 'pako'

const urlExample = document.getElementById('link-example').href
var exampleData
var view = 'samples'

$('#mainTab a').on('click', function(e) {
  e.preventDefault()
  $(this).tab('show')
})

const fileInput = document.getElementById('inputFile')
const resultLink = document.getElementById('link-results')
const resultInfo = document.getElementById('result-info')
const sampleSelect = document.getElementById('selectSample')
const resultContainer = document.getElementById('result-container')
const chartsContainer = document.getElementById('charts-container')

const submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', function() {
  resultLink.click()
  run()
})

const exampleButton = document.getElementById('btn-example')
exampleButton.addEventListener('click', showExample)

const trellisSampleButton = document.getElementById('btn-trellis-sample')
trellisSampleButton.addEventListener('click', trellisSample)

const trellisChromosomeButton = document.getElementById(
  'btn-trellis-chromosome'
)
trellisChromosomeButton.addEventListener('click', trellisChromosome)

sampleSelect.addEventListener('change', handleChangeSample)

var inputData, samples, chromosomes

function run() {
  const file = fileInput.files[0]

  // TODO error handling
  if (file === undefined) return

  view = 'samples'
  trellisSampleButton.classList.add('active')
  trellisChromosomeButton.classList.remove('active')
  hideElement(resultContainer)
  chartsContainer.innerHTML = ''
  showElement(resultInfo)

  const fileReader = new FileReader()
  const isGzip = /\.gz$/.test(file.name)
  if (isGzip) {
    fileReader.readAsArrayBuffer(file)
  } else {
    fileReader.readAsText(file)
  }
  fileReader.onload = event => {
    let content = event.target.result
    if (isGzip) {
      content = pako.ungzip(content, { to: 'string' })
    }
    inputData = JSON.parse(content).data
    setupSamples()
    setupChromosomes()
    setupSelect()
    setTimeout(() => {
      vis(0)
    }, 0)
  }
}

function setupSamples() {
  samples = inputData.map(rec => rec.sample)
}

function setupChromosomes() {
  chromosomes = Array.from(
    new Set(
      inputData
        .map(rec => rec.coverages.map(cov => cov.chromosome))
        .reduce((chroms, chrom) => chroms.concat(chrom))
    )
  )
}

function setupSelect() {
  if (view === 'samples') {
    sampleSelect.innerHTML = samples
      .map((sample, index) => `<option value="${index}">${sample}</option>`)
      .join('\n')
  } else {
    sampleSelect.innerHTML = chromosomes
      .map(
        (chromosome, index) => `<option value="${index}">${chromosome}</option>`
      )
      .join('\n')
  }
}

function trellisSample() {
  if (view === 'samples') {
    return
  }
  trellisSampleButton.classList.add('active')
  trellisChromosomeButton.classList.remove('active')
  view = 'samples'
  hideElement(resultContainer)
  chartsContainer.innerHTML = ''
  showElement(resultInfo)
  setupSelect()
  setTimeout(() => {
    vis(0)
  }, 0)
}

function trellisChromosome() {
  if (view === 'chromosomes') {
    return
  }
  trellisChromosomeButton.classList.add('active')
  trellisSampleButton.classList.remove('active')
  view = 'chromosomes'
  hideElement(resultContainer)
  chartsContainer.innerHTML = ''
  showElement(resultInfo)
  setupSelect()
  setTimeout(() => {
    vis(0)
  }, 0)
}

function handleChangeSample(event) {
  const index = Number.parseInt(event.target.value)
  hideElement(resultContainer)
  chartsContainer.innerHTML = ''
  showElement(resultInfo)
  setTimeout(() => {
    vis(index)
  }, 0)
}

function getData(index) {
  const ret = []
  if (view === 'samples') {
    const sample = samples[index]
    const slot = inputData.find(rec => rec.sample === sample)
    for (const coverage of slot.coverages) {
      ret.push({
        sample,
        chromosome: coverage.chromosome,
        counts: coverage.counts,
        positions: coverage.positions
      })
    }
  } else {
    const chromosome = chromosomes[index]
    for (const rec of inputData) {
      for (const coverage of rec.coverages.filter(
        coverage => coverage.chromosome === chromosome
      )) {
        ret.push({
          sample: rec.sample,
          chromosome,
          counts: coverage.counts,
          positions: coverage.positions
        })
      }
    }
  }
  return ret
}

function vis(index) {
  const data = getData(index)
  showElement(resultContainer)
  for (const rec of data) {
    const coverageContainer = document.createElement('div')
    const chartContainer = document.createElement('div')
    coverageContainer.appendChild(chartContainer)
    chartsContainer.appendChild(coverageContainer)
    const title = `${rec.chromosome} – ${rec.sample}`
    renderStrandedCoverageChart(
      chartContainer,
      { positions: rec.positions, counts: rec.counts },
      title
    )
  }
  hideElement(resultInfo)
}

function renderStrandedCoverageChart(container, data, title) {
  const trace1 = {
    x: data.positions,
    y: data.counts[0].values,
    name: data.counts[0].label,
    type: 'bar'
  }

  const trace2 = {
    x: data.positions,
    y: data.counts[1].values,
    yaxis: 'y2',
    name: data.counts[1].label,
    type: 'bar'
  }

  const traces = [trace1, trace2]

  const layout = {
    title: title || '',
    yaxis: {
      title: 'count',
      domain: [0.54, 1]
    },
    yaxis2: {
      domain: [0, 0.46],
      autorange: 'reversed'
    }
  }

  Plotly.newPlot(container, traces, layout)
}

function showExample() {
  trellisSampleButton.classList.add('active')
  trellisChromosomeButton.classList.remove('active')
  view = 'samples'

  hideElement(resultContainer)
  chartsContainer.innerHTML = ''
  showElement(resultInfo)
  resultLink.click()

  if (exampleData) {
    inputData = exampleData
    setupSamples()
    setupChromosomes()
    setupSelect()
    setTimeout(() => {
      vis(0)
    }, 0)
  } else {
    axios
      .get(urlExample, {
        responseType: 'arraybuffer'
      })
      .then(response => {
        const content = pako.ungzip(response.data, { to: 'string' })
        exampleData = JSON.parse(content).data
        inputData = exampleData
        setupSamples()
        setupChromosomes()
        setupSelect()
        setTimeout(() => {
          vis(0)
        }, 0)
      })
      .catch(error => {
        // FIXME proper error handling
        console.error(error)
      })
  }
}

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}
